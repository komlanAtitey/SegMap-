#!/usr/bin/env python3
"""
synthetic_tissue.py
===================

A ground-truth generative model of spatial transcriptomic tissue for
benchmarking segmentation, transcript assignment, marker calling and
denoising. Implements the five-stage design:

  1. Cell-type transcriptional programs and exclusive marker panels.
  2. Type-specific cell and nucleus morphology (size, shape, nuclear position).
  3. Subcellular transcript localization (nucleus / cytosol / whole-cell).
  4. Non-overlapping packing into a confined tissue domain (with optional,
     controlled overlap and interstitial gaps).
  5. Low-level background (noise) transcripts belonging to no cell.

Everything that is unobservable in a real experiment is retained as ground
truth: each molecule's true cell of origin (-1 = background), each cell's
type, and each cell/nucleus polygon. Outputs are written in a transcript-table
format compatible with downstream pipelines, plus the cell table, gene panel,
polygon archive, a rasterized label mask, and an overview figure.

Dependencies: numpy, pandas, matplotlib.

Example
-------
python synthetic_tissue.py --n-cells 150 --n-genes 200 --n-types 6 \
    --out-dir synthetic_out --seed 0
"""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================
# Configuration
# ============================================================

@dataclass
class TissueConfig:
    n_cells: int = 150
    n_genes: int = 200
    n_types: int = 6                 # 6-7 recommended
    markers_per_type: int = 12       # 10-15 exclusive markers per type
    expr_genes_min: int = 50         # expressed genes per type (incl. markers)
    expr_genes_max: int = 100

    # expression means (expected counts per cell)
    mu_marker: Tuple[float, float] = (8.0, 18.0)     # high, exclusive
    mu_shared: Tuple[float, float] = (0.5, 3.0)      # baseline expression level
    leak: float = 0.0                                # mean for non-expressed genes
    # marker_mode: "distinct" = non-marker genes share one global baseline across
    #   ALL types, so exclusive markers are the ONLY differentially expressed genes
    #   (crisp marker ground truth; marker precision/F1 are meaningful).
    # "graded" = each type gets its own random means for a random shared subset
    #   (more realistic non-marker variation, but marker F1 is capped because many
    #   shared genes are legitimately DE and are not in the exclusive marker set).
    marker_mode: str = "distinct"

    # subcellular compartment mix over the panel
    frac_nucleus: float = 0.25
    frac_cytosol: float = 0.45
    # remainder -> whole-cell

    # count model
    count_model: str = "nb"          # "poisson" or "nb"
    nb_dispersion: float = 5.0        # NB shape r; var = mu + mu^2/r

    # packing
    areal_density: float = 0.55       # target fraction of domain covered by cells
    spacing_factor: float = 1.0       # >=1 abut with gaps; <1 forces contact
    overlap_frac: float = 0.08        # fraction of cells allowed to overlap a neighbour
    overlap_spacing: float = 0.70     # relaxed spacing for the overlapping fraction
    max_place_attempts: int = 2000

    # background noise
    lambda_bg: float = 0.03           # background transcripts per unit area

    # geometry discretization
    n_boundary: int = 48              # polygon vertices per cell
    nucleus_scale: Tuple[float, float] = (0.40, 0.55)  # nucleus size vs cell

    seed: int = 0


# Type-specific morphology templates (semi-axis ranges in arbitrary length
# units, irregularity amplitude, and nuclear eccentricity 0=central..1=polar).
def default_type_templates(n_types: int, rng: np.random.Generator) -> List[dict]:
    base = [
        dict(name="A", a=(18, 24), b=(7, 10),  irregularity=0.10, nuc_ecc=0.15),  # large, elongated
        dict(name="B", a=(10, 14), b=(9, 13),  irregularity=0.28, nuc_ecc=0.55),  # intermediate, irregular
        dict(name="C", a=(6, 9),   b=(5, 7),   irregularity=0.08, nuc_ecc=0.10),  # small, ellipsoidal
        dict(name="D", a=(10, 13), b=(9, 12),  irregularity=0.12, nuc_ecc=0.20),  # medium, round
        dict(name="E", a=(11, 15), b=(4, 6),   irregularity=0.15, nuc_ecc=0.60),  # small, elongated
        dict(name="F", a=(15, 19), b=(13, 17), irregularity=0.22, nuc_ecc=0.25),  # large, round-irregular
        dict(name="G", a=(5, 7),   b=(4, 6),   irregularity=0.10, nuc_ecc=0.15),  # tiny
    ]
    return base[:n_types]


# ============================================================
# Geometry helpers (numpy-only; no shapely)
# ============================================================

def _ellipse_radius(phi: np.ndarray, a: float, b: float) -> np.ndarray:
    """Radius of an axis-aligned ellipse (semi-axes a, b) at polar angle phi."""
    return (a * b) / np.sqrt((b * np.cos(phi)) ** 2 + (a * np.sin(phi)) ** 2)


def make_polygon(cx, cy, a, b, theta, irregularity, n_vert, rng,
                 n_harmonics=4) -> np.ndarray:
    """Star-shaped polygon: an ellipse whose radius is perturbed by a few
    low-order harmonics. Radially monotone, so point-in-polygon and rejection
    sampling are well behaved."""
    phi = np.linspace(0, 2 * np.pi, n_vert, endpoint=False)
    perturb = np.ones_like(phi)
    for m in range(2, 2 + n_harmonics):
        amp = irregularity * rng.uniform(0.3, 1.0) / (m - 1)
        perturb += amp * np.cos(m * phi + rng.uniform(0, 2 * np.pi))
    perturb = np.clip(perturb, 0.35, None)
    r = _ellipse_radius(phi, a, b) * perturb
    xl, yl = r * np.cos(phi), r * np.sin(phi)
    ct, st = np.cos(theta), np.sin(theta)
    x = cx + ct * xl - st * yl
    y = cy + st * xl + ct * yl
    return np.column_stack([x, y])


def points_in_poly(pts: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Vectorized ray-casting point-in-polygon for many points, one polygon."""
    x, y = pts[:, 0], pts[:, 1]
    vx, vy = poly[:, 0], poly[:, 1]
    n = len(poly)
    inside = np.zeros(len(pts), dtype=bool)
    j = n - 1
    for i in range(n):
        xi, yi, xj, yj = vx[i], vy[i], vx[j], vy[j]
        cond = ((yi > y) != (yj > y)) & (
            x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi
        )
        inside ^= cond
        j = i
    return inside


def sample_in_region(poly: np.ndarray, n: int, rng: np.random.Generator,
                     exclude: Optional[np.ndarray] = None,
                     max_iter: int = 200) -> np.ndarray:
    """Rejection-sample n points inside `poly` (optionally excluding `exclude`)."""
    if n <= 0:
        return np.empty((0, 2))
    mn = poly.min(0)
    mx = poly.max(0)
    out = np.empty((0, 2))
    it = 0
    while len(out) < n and it < max_iter:
        batch = max(2 * (n - len(out)), 64)
        cand = rng.uniform(mn, mx, size=(batch, 2))
        ok = points_in_poly(cand, poly)
        if exclude is not None:
            ok &= ~points_in_poly(cand, exclude)
        if ok.any():
            out = np.vstack([out, cand[ok]])
        it += 1
    return out[:n]


# ============================================================
# Gene panel and expression program
# ============================================================

def build_panel(cfg: TissueConfig, rng: np.random.Generator):
    """Assign exclusive markers, per-type expressed sets, compartments, and the
    type x gene expression-mean matrix MU."""
    G, T = cfg.n_genes, cfg.n_types
    m = cfg.markers_per_type
    total_markers = m * T
    if total_markers > G:
        raise ValueError(f"Need {total_markers} marker genes but panel has {G}.")

    perm = rng.permutation(G)
    marker_ids = perm[:total_markers]
    type_markers = {t: marker_ids[t * m:(t + 1) * m].tolist() for t in range(T)}
    shared_pool = perm[total_markers:].tolist()   # non-marker genes (shared)

    # compartment per gene
    comp = np.empty(G, dtype=object)
    r = rng.random(G)
    comp[r < cfg.frac_nucleus] = "nucleus"
    comp[(r >= cfg.frac_nucleus) & (r < cfg.frac_nucleus + cfg.frac_cytosol)] = "cytosol"
    comp[r >= cfg.frac_nucleus + cfg.frac_cytosol] = "whole"

    MU = np.full((T, G), cfg.leak, dtype=float)
    expressed: Dict[int, List[int]] = {}

    if cfg.marker_mode == "distinct":
        # One global baseline per gene; a single shared set expressed EQUALLY by
        # all types, so only the exclusive markers are differentially expressed.
        base = rng.uniform(cfg.mu_shared[0], cfg.mu_shared[1], size=G)
        target = int(round(0.5 * (cfg.expr_genes_min + cfg.expr_genes_max)))
        n_shared = int(np.clip(target - m, 1, len(shared_pool)))
        common_shared = np.asarray(shared_pool[:n_shared], dtype=int)
        for t in range(T):
            mk = type_markers[t]
            MU[t, common_shared] = base[common_shared]          # identical across types
            MU[t, mk] = rng.uniform(*cfg.mu_marker, size=len(mk))  # exclusive, high
            expressed[t] = mk + common_shared.tolist()
    elif cfg.marker_mode == "graded":
        for t in range(T):
            mk = type_markers[t]
            MU[t, mk] = rng.uniform(*cfg.mu_marker, size=len(mk))
            n_expr = int(rng.integers(cfg.expr_genes_min, cfg.expr_genes_max + 1))
            n_shared = max(0, n_expr - len(mk))
            shared_sel = rng.choice(shared_pool, size=min(n_shared, len(shared_pool)),
                                    replace=False)
            MU[t, shared_sel] = rng.uniform(*cfg.mu_shared, size=len(shared_sel))
            expressed[t] = mk + shared_sel.tolist()
    else:
        raise ValueError(f"Unknown marker_mode: {cfg.marker_mode}")

    gene_names = [f"Gene_{g:04d}" for g in range(G)]
    panel = pd.DataFrame({
        "gene": gene_names,
        "gene_index": np.arange(G),
        "compartment": comp,
        "is_marker": [any(g in type_markers[t] for t in range(T)) for g in range(G)],
        "marker_of_type": [next((t for t in range(T) if g in type_markers[t]), -1)
                           for g in range(G)],
    })
    return MU, comp, gene_names, type_markers, expressed, panel


# ============================================================
# Cell specification and packing
# ============================================================

def draw_cell_specs(cfg: TissueConfig, templates, rng, mix=None):
    """Draw per-cell type and size BEFORE placement, so the domain can be sized."""
    T = cfg.n_types
    if mix is None:
        mix = np.ones(T) / T
    types = rng.choice(T, size=cfg.n_cells, p=mix)
    specs = []
    for t in types:
        tpl = templates[t]
        a = rng.uniform(*tpl["a"])
        b = rng.uniform(*tpl["b"])
        theta = rng.uniform(0, np.pi)
        irow = tpl["irregularity"]
        r_eff = float(np.sqrt(a * b))
        r_bound = float(max(a, b) * (1 + irow))
        specs.append(dict(type=int(t), a=a, b=b, theta=theta,
                          irregularity=irow, nuc_ecc=tpl["nuc_ecc"],
                          r_eff=r_eff, r_bound=r_bound))
    return specs


def pack_cells(specs, cfg: TissueConfig, rng) -> Tuple[np.ndarray, List[int]]:
    """Dart-throwing packing. Returns placed centroids and the kept spec indices."""
    mean_area = np.mean([np.pi * s["a"] * s["b"] for s in specs])
    domain_area = cfg.n_cells * mean_area / cfg.areal_density
    W = H = float(np.sqrt(domain_area))

    order = np.argsort([-s["r_eff"] for s in specs])  # largest first
    centers = np.full((len(specs), 2), np.nan)
    placed: List[int] = []
    for idx in order:
        s = specs[idx]
        rb = s["r_bound"]
        allow_overlap = rng.random() < cfg.overlap_frac
        sf = cfg.overlap_spacing if allow_overlap else cfg.spacing_factor
        for _ in range(cfg.max_place_attempts):
            c = rng.uniform([rb, rb], [W - rb, H - rb])
            if not placed:
                centers[idx] = c
                placed.append(idx)
                break
            P = centers[placed]
            d = np.hypot(P[:, 0] - c[0], P[:, 1] - c[1])
            need = sf * (np.array([specs[j]["r_eff"] for j in placed]) + s["r_eff"])
            if np.all(d >= need):
                centers[idx] = c
                placed.append(idx)
                break
    placed_sorted = sorted(placed)
    return centers, placed_sorted, (W, H)


# ============================================================
# Transcript generation
# ============================================================

def _counts(mu_row: np.ndarray, cfg: TissueConfig, rng) -> np.ndarray:
    if cfg.count_model == "poisson":
        return rng.poisson(mu_row)
    r = cfg.nb_dispersion
    mu = np.maximum(mu_row, 1e-9)
    p = r / (r + mu)
    out = np.zeros_like(mu, dtype=int)
    nz = mu_row > 0
    out[nz] = rng.negative_binomial(r, p[nz])
    return out


def generate_transcripts(cells: pd.DataFrame, polys, nuclei, MU, comp,
                         gene_names, cfg: TissueConfig, rng):
    """Emit one row per transcript with its true cell_id, type, gene, compartment."""
    comp_arr = np.asarray(comp)
    recs_x, recs_y, recs_g, recs_c, recs_t, recs_comp = [], [], [], [], [], []
    for _, row in cells.iterrows():
        cid = int(row["cell_id"])
        t = int(row["cell_type"])
        cell_poly = polys[cid]
        nuc_poly = nuclei[cid]
        counts = _counts(MU[t], cfg, rng)
        if counts.sum() == 0:
            continue
        for region_name, region_poly, excl in (
            ("nucleus", nuc_poly, None),
            ("cytosol", cell_poly, nuc_poly),
            ("whole", cell_poly, None),
        ):
            gsel = np.where((comp_arr == region_name) & (counts > 0))[0]
            if len(gsel) == 0:
                continue
            gene_rep = np.repeat(gsel, counts[gsel])
            pts = sample_in_region(region_poly, len(gene_rep), rng, exclude=excl)
            if len(pts) < len(gene_rep):          # rare: tight exclusion
                gene_rep = gene_rep[:len(pts)]
            recs_x.append(pts[:, 0]); recs_y.append(pts[:, 1])
            recs_g.append(gene_rep)
            recs_c.append(np.full(len(gene_rep), cid))
            recs_t.append(np.full(len(gene_rep), t))
            recs_comp.append(np.full(len(gene_rep), region_name, dtype=object))

    if not recs_x:
        return pd.DataFrame(columns=["x", "y", "gene", "cell_id", "cell_type", "compartment"])
    df = pd.DataFrame({
        "x": np.concatenate(recs_x),
        "y": np.concatenate(recs_y),
        "gene_index": np.concatenate(recs_g),
        "cell_id": np.concatenate(recs_c),
        "cell_type": np.concatenate(recs_t),
        "compartment": np.concatenate(recs_comp),
    })
    df["gene"] = [gene_names[g] for g in df["gene_index"]]
    return df


def generate_background(domain, n_genes, gene_names, cfg: TissueConfig, rng,
                        gene_weights=None):
    W, H = domain
    n_bg = rng.poisson(cfg.lambda_bg * W * H)
    if n_bg == 0:
        return pd.DataFrame(columns=["x", "y", "gene", "cell_id", "cell_type",
                                     "compartment", "gene_index"])
    x = rng.uniform(0, W, n_bg)
    y = rng.uniform(0, H, n_bg)
    g = rng.choice(n_genes, size=n_bg, p=gene_weights)
    return pd.DataFrame({
        "x": x, "y": y, "gene_index": g,
        "gene": [gene_names[i] for i in g],
        "cell_id": -1, "cell_type": -1, "compartment": "background",
    })


# ============================================================
# Rasterized ground-truth label mask
# ============================================================

def rasterize_labels(polys, nuclei, domain, pixel_size, cell_ids):
    W, H = domain
    nx, ny = int(np.ceil(W / pixel_size)), int(np.ceil(H / pixel_size))
    cell_lab = np.zeros((ny, nx), dtype=np.int32)      # 0 = background
    nuc_lab = np.zeros((ny, nx), dtype=np.int32)
    for cid in cell_ids:
        for poly, lab in ((polys[cid], cell_lab), (nuclei[cid], nuc_lab)):
            mnx, mny = poly.min(0); mxx, mxy = poly.max(0)
            i0, i1 = max(0, int(mnx / pixel_size)), min(nx, int(np.ceil(mxx / pixel_size)))
            j0, j1 = max(0, int(mny / pixel_size)), min(ny, int(np.ceil(mxy / pixel_size)))
            if i1 <= i0 or j1 <= j0:
                continue
            xs = (np.arange(i0, i1) + 0.5) * pixel_size
            ys = (np.arange(j0, j1) + 0.5) * pixel_size
            gx, gy = np.meshgrid(xs, ys)
            pts = np.column_stack([gx.ravel(), gy.ravel()])
            inside = points_in_poly(pts, poly).reshape(gy.shape)
            sub = lab[j0:j1, i0:i1]
            sub[inside] = cid + 1                      # labels are 1-based
    return cell_lab, nuc_lab


# ============================================================
# Build + save + verify
# ============================================================

def build_tissue(cfg: TissueConfig):
    rng = np.random.default_rng(cfg.seed)
    templates = default_type_templates(cfg.n_types, rng)
    MU, comp, gene_names, type_markers, expressed, panel = build_panel(cfg, rng)

    specs = draw_cell_specs(cfg, templates, rng)
    centers, kept, domain = pack_cells(specs, cfg, rng)

    polys, nuclei, rows = {}, {}, []
    for new_id, idx in enumerate(kept):
        s = specs[idx]
        cx, cy = centers[idx]
        poly = make_polygon(cx, cy, s["a"], s["b"], s["theta"],
                            s["irregularity"], cfg.n_boundary, rng)
        # nucleus: smaller, rounder, offset toward a random pole by nuc_ecc
        nscale = rng.uniform(*cfg.nucleus_scale)
        na, nb = s["a"] * nscale, s["b"] * nscale
        maxoff = (1 - nscale) * min(s["a"], s["b"]) * 0.8
        ang = rng.uniform(0, 2 * np.pi)
        off = s["nuc_ecc"] * maxoff
        ncx = cx + off * np.cos(ang + s["theta"])
        ncy = cy + off * np.sin(ang + s["theta"])
        nuc = make_polygon(ncx, ncy, na, nb, s["theta"],
                           s["irregularity"] * 0.5, cfg.n_boundary, rng)
        polys[new_id] = poly
        nuclei[new_id] = nuc
        rows.append(dict(cell_id=new_id, cell_type=s["type"],
                         type_name=templates[s["type"]]["name"],
                         cx=cx, cy=cy, nucleus_cx=ncx, nucleus_cy=ncy,
                         semi_major=s["a"], semi_minor=s["b"],
                         orientation=s["theta"], irregularity=s["irregularity"],
                         area_approx=float(np.pi * s["a"] * s["b"])))
    cells = pd.DataFrame(rows)

    tx = generate_transcripts(cells, polys, nuclei, MU, comp, gene_names, cfg, rng)
    bg = generate_background(domain, cfg.n_genes, gene_names, cfg, rng)
    transcripts = pd.concat([tx, bg], ignore_index=True)
    transcripts = transcripts[["x", "y", "gene", "gene_index", "cell_id",
                               "cell_type", "compartment"]]

    return dict(cfg=cfg, domain=domain, cells=cells, transcripts=transcripts,
                polys=polys, nuclei=nuclei, panel=panel, MU=MU,
                type_markers=type_markers, templates=templates)


def verify(result) -> dict:
    """Ground-truth sanity checks."""
    tx = result["transcripts"]
    polys = result["polys"]
    real = tx[tx["cell_id"] >= 0]
    # fraction of real transcripts inside their own assigned cell polygon
    inside = 0
    for cid, grp in real.groupby("cell_id"):
        pts = grp[["x", "y"]].to_numpy()
        inside += int(points_in_poly(pts, polys[int(cid)]).sum())
    frac_inside = inside / max(len(real), 1)
    # pairwise polygon overlap estimate via bounding circles actually touching
    cells = result["cells"]
    cxy = cells[["cx", "cy"]].to_numpy()
    reff = np.sqrt(cells["semi_major"] * cells["semi_minor"]).to_numpy()
    n = len(cells); overlaps = 0; pairs = 0
    for i in range(n):
        d = np.hypot(cxy[i + 1:, 0] - cxy[i, 0], cxy[i + 1:, 1] - cxy[i, 1])
        overlaps += int(np.sum(d < (reff[i + 1:] + reff[i]) * 0.85))
        pairs += len(d)
    return dict(
        n_cells=int(len(cells)),
        n_transcripts=int(len(tx)),
        n_real=int(len(real)),
        n_background=int((tx["cell_id"] < 0).sum()),
        background_fraction=float((tx["cell_id"] < 0).mean()),
        frac_real_inside_true_cell=round(frac_inside, 6),
        approx_overlapping_pairs=overlaps,
        total_pairs=pairs,
    )


def save_outputs(result, out_dir, pixel_size=0.5, make_figure=True):
    os.makedirs(out_dir, exist_ok=True)
    cfg = result["cfg"]
    result["transcripts"].to_csv(os.path.join(out_dir, "transcripts.csv"), index=False)
    result["cells"].to_csv(os.path.join(out_dir, "cells.csv"), index=False)
    result["panel"].to_csv(os.path.join(out_dir, "gene_panel.csv"), index=False)
    np.save(os.path.join(out_dir, "expression_mu.npy"), result["MU"])

    # polygons
    np.savez_compressed(
        os.path.join(out_dir, "polygons.npz"),
        **{f"cell_{k}": v for k, v in result["polys"].items()},
        **{f"nucleus_{k}": v for k, v in result["nuclei"].items()},
    )
    # label mask
    cell_lab, nuc_lab = rasterize_labels(
        result["polys"], result["nuclei"], result["domain"], pixel_size,
        list(result["polys"].keys()))
    np.save(os.path.join(out_dir, "label_mask_cells.npy"), cell_lab)
    np.save(os.path.join(out_dir, "label_mask_nuclei.npy"), nuc_lab)

    stats = verify(result)
    stats["domain"] = list(result["domain"])
    stats["pixel_size"] = pixel_size
    with open(os.path.join(out_dir, "ground_truth_summary.json"), "w") as f:
        json.dump({"config": asdict(cfg), "stats": stats}, f, indent=2, default=float)

    if make_figure:
        _overview_figure(result, os.path.join(out_dir, "overview.png"))
    return stats


def _select_zoom_window(cells, target_cells=9):
    """Choose a square window centred on a cell that maximizes cell-type
    diversity (then cell count) within it, so the zoom shows a representative,
    contiguous neighbourhood rather than scattered cells."""
    cxy = cells[["cx", "cy"]].to_numpy()
    types = cells["cell_type"].to_numpy()
    mean_area = float(cells["area_approx"].mean())
    side = float(np.sqrt(target_cells * mean_area / 0.5))   # ~target_cells cells
    half = side / 2.0
    best = None
    for i in range(len(cells)):
        cx, cy = cxy[i]
        m = (np.abs(cxy[:, 0] - cx) <= half) & (np.abs(cxy[:, 1] - cy) <= half)
        n = int(m.sum()); ndist = int(len(np.unique(types[m])))
        if n < 4:
            continue
        # prefer diversity, then a count close to target (not too crowded)
        score = ndist * 1000 - abs(n - target_cells)
        if best is None or score > best[0]:
            best = (score, cx, cy)
    if best is None:
        cx, cy = cxy.mean(0)
    else:
        _, cx, cy = best
    return cx, cy, half


def _overview_figure(result, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPoly, Rectangle, ConnectionPatch
    from matplotlib.collections import PatchCollection

    cfg = result["cfg"]; W, H = result["domain"]
    cells = result["cells"]; polys = result["polys"]; nuclei = result["nuclei"]
    tx = result["transcripts"]
    cmap = plt.get_cmap("tab10")
    comp_color = {"nucleus": "#d62728", "cytosol": "#1f77b4", "whole": "#2ca02c"}

    # choose the zoom window and the cells whose centroids fall inside it
    zx, zy, half = _select_zoom_window(cells)
    x0, x1, y0, y1 = zx - half, zx + half, zy - half, zy + half
    in_win = cells[(cells.cx >= x0) & (cells.cx <= x1) &
                   (cells.cy >= y0) & (cells.cy <= y1)]
    win_ids = set(in_win.cell_id.astype(int))
    n_types_win = in_win.cell_type.nunique()

    fig, ax = plt.subplots(1, 2, figsize=(15.5, 7.4))
    axL, axR = ax

    # ---- (a) whole tissue coloured by cell type ----
    patches, colors = [], []
    for _, r in cells.iterrows():
        patches.append(MplPoly(polys[int(r.cell_id)], closed=True))
        colors.append(cmap(int(r.cell_type) % 10))
    axL.add_collection(PatchCollection(patches, facecolor=colors, edgecolor="k",
                                       linewidths=0.3, alpha=0.75))
    for _, r in cells.iterrows():
        axL.add_patch(MplPoly(nuclei[int(r.cell_id)], closed=True,
                              facecolor="white", edgecolor="0.3",
                              linewidth=0.3, alpha=0.6))
    bg = tx[tx.cell_id < 0]
    axL.scatter(bg.x, bg.y, s=1, c="0.5", alpha=0.4)
    # zoom rectangle
    axL.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                            edgecolor="k", linewidth=1.8, zorder=5))
    axL.set_xlim(0, W); axL.set_ylim(0, H); axL.set_aspect("equal")
    axL.set_title(f"Synthetic tissue: {len(cells)} cells, {cfg.n_types} types\n"
                  f"nuclei in white; grey = background noise; box = zoom region",
                  fontsize=11, loc="left")
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=cmap(t % 10), markersize=8,
                          label=result["templates"][t]["name"])
               for t in range(cfg.n_types)]
    axL.legend(handles=handles, title="cell type", loc="upper right",
               fontsize=8, framealpha=0.9)

    # ---- (b) zoom: same cells, light type fill + compartment-coloured points ----
    for _, r in in_win.iterrows():
        cid = int(r.cell_id)
        axR.add_patch(MplPoly(polys[cid], closed=True,
                              facecolor=cmap(int(r.cell_type) % 10),
                              edgecolor=cmap(int(r.cell_type) % 10),
                              linewidth=1.0, alpha=0.18))
        axR.add_patch(MplPoly(nuclei[cid], closed=True, facecolor="none",
                              edgecolor="0.4", linewidth=0.9, linestyle="--"))
    sel = tx[tx.cell_id.isin(win_ids)]
    for comp_name, col in comp_color.items():
        s = sel[sel.compartment == comp_name]
        axR.scatter(s.x, s.y, s=6, c=col, alpha=0.8, label=comp_name,
                    edgecolors="none")
    # background noise inside the window too
    bgw = bg[(bg.x >= x0) & (bg.x <= x1) & (bg.y >= y0) & (bg.y <= y1)]
    axR.scatter(bgw.x, bgw.y, s=5, c="0.6", alpha=0.5, marker="x",
                label="background")
    axR.set_xlim(x0, x1); axR.set_ylim(y0, y1); axR.set_aspect("equal")
    axR.set_title(f"Subcellular transcript localization (zoom)\n"
                  f"{len(in_win)} cells, {n_types_win} types; "
                  f"nucleus outline dashed", fontsize=11, loc="left")
    axR.legend(title="compartment", loc="upper right", fontsize=8, framealpha=0.9)

    # ---- connector lines linking the rectangle to the zoom panel ----
    for (xa, ya), (xb, yb) in [((x1, y1), (x0, y1)), ((x1, y0), (x0, y0))]:
        con = ConnectionPatch(xyA=(xa, ya), coordsA=axL.transData,
                              xyB=(xb, yb), coordsB=axR.transData,
                              color="0.3", linewidth=1.0, linestyle="-",
                              zorder=10)
        fig.add_artist(con)

    plt.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# CLI
# ============================================================

def build_parser():
    p = argparse.ArgumentParser(description="Generate a ground-truth synthetic "
                                            "spatial-transcriptomics tissue.")
    p.add_argument("--n-cells", type=int, default=150)
    p.add_argument("--n-genes", type=int, default=200)
    p.add_argument("--n-types", type=int, default=6)
    p.add_argument("--markers-per-type", type=int, default=12)
    p.add_argument("--marker-mode", choices=["distinct", "graded"], default="distinct",
                   help="distinct = crisp exclusive markers (recommended for "
                        "marker-recovery scoring); graded = realistic shared-gene "
                        "variation (marker F1 capped).")
    p.add_argument("--count-model", choices=["poisson", "nb"], default="nb")
    p.add_argument("--nb-dispersion", type=float, default=5.0)
    p.add_argument("--areal-density", type=float, default=0.55)
    p.add_argument("--overlap-frac", type=float, default=0.08)
    p.add_argument("--lambda-bg", type=float, default=0.03)
    p.add_argument("--pixel-size", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out-dir", default="synthetic_out")
    p.add_argument("--no-figure", action="store_true")
    return p


def main():
    a = build_parser().parse_args()
    cfg = TissueConfig(
        n_cells=a.n_cells, n_genes=a.n_genes, n_types=a.n_types,
        markers_per_type=a.markers_per_type, marker_mode=a.marker_mode,
        count_model=a.count_model,
        nb_dispersion=a.nb_dispersion, areal_density=a.areal_density,
        overlap_frac=a.overlap_frac, lambda_bg=a.lambda_bg, seed=a.seed,
    )
    print("Building synthetic tissue...")
    result = build_tissue(cfg)
    stats = save_outputs(result, a.out_dir, pixel_size=a.pixel_size,
                         make_figure=not a.no_figure)
    print("\nGround-truth summary:")
    for k, v in stats.items():
        print(f"  {k}: {v}")
    print(f"\nOutputs written to: {os.path.abspath(a.out_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
