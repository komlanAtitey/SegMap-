"""
refine_boundaries.py — SegMap transcriptome-guided cell-boundary refinement (Fig. 1a).

Option 1: transcript-level reassignment at cell borders, initialized from the
nucleus-anchored segmentation. Only boundary-associated transcripts are
reconsidered; nuclear architecture and cell morphology are preserved. This is an
optional post-assignment step — it refines the transcript -> cell mapping and the
cell-by-gene matrix, and it does NOT move pixel masks, re-segment nuclei, or alter
cell centroids.

Boundary zone. A transcript is boundary-associated if its distances to the nearest
and second-nearest nuclei satisfy (d2 - d1) / d1 < tau  (tau = ambiguity threshold).

Scoring. For a boundary transcript t of gene g and candidate cell i:
    l(t -> i) = log theta_{i,g} - lambda * d(t, n_i)
with theta_i the (confidently-assigned) expression profile of cell i and
d(t, n_i) the distance from t to cell i's nucleus centroid. Assignment is argmax.

EM loop. (1) estimate profiles from confidently-assigned transcripts;
(2) reassign boundary transcripts by max score; (3) recompute; repeat to
convergence or max_iter.

Validation (returned + written). Primary metric: multi-lineage contamination
    C_i = 1 - max_k p_{ik},   C = mean_i C_i,   Delta C = C_before - C_after (>0 = better),
where p_{ik} is the fraction of cell i's transcripts from lineage k, given a
user-supplied gene->lineage map. Refinement should be claimed only if Delta C > 0.

Off by default in the pipeline; writes to a separate folder; never overwrites the
original assignment or cell-by-gene matrix.

Dependencies: numpy, pandas, scipy.
"""
from __future__ import annotations
import os, json
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

__all__ = ["refine_boundaries", "contamination_score", "validate_refinement", "derive_gene_lineage_from_clusters", "run_refinement"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _profiles(counts: csr_matrix, cell_index: np.ndarray, n_cells: int,
              pseudocount: float = 1.0) -> np.ndarray:
    """Per-cell CP-normalized expression profile theta_i (n_cells x G) from the
    transcript->cell assignment `cell_index` (row per transcript, value = cell row
    or -1). Uses a Dirichlet-style pseudocount so log theta is finite."""
    G = counts.shape[1]
    M = np.full((n_cells, G), 0.0)
    valid = cell_index >= 0
    ci = cell_index[valid]
    # counts is transcript x gene one-hot; accumulate per cell
    np.add.at(M, ci, counts[valid].toarray())
    M += pseudocount
    M /= M.sum(1, keepdims=True)
    return M


def _nearest_two_nuclei(tx_xy: np.ndarray, centroids_xy: np.ndarray):
    """Return (i1, d1, i2, d2): indices/distances of the nearest and second-nearest
    nucleus centroid for every transcript."""
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=2).fit(centroids_xy)
    dist, idx = nn.kneighbors(tx_xy)
    return idx[:, 0], dist[:, 0], idx[:, 1], dist[:, 1]


def contamination_score(cell_index: np.ndarray, gene_index: np.ndarray,
                        gene_lineage: np.ndarray, n_cells: int):
    """Multi-lineage contamination. gene_lineage maps each gene -> lineage id
    (or -1 to ignore). Returns (C_mean, C_per_cell, n_cells_scored)."""
    lut = gene_lineage[gene_index]                       # lineage per transcript
    valid = (cell_index >= 0) & (lut >= 0)
    ci, li = cell_index[valid], lut[valid]
    K = int(gene_lineage[gene_lineage >= 0].max()) + 1 if (gene_lineage >= 0).any() else 0
    if K == 0:
        return np.nan, np.full(n_cells, np.nan), 0
    counts = np.zeros((n_cells, K))
    np.add.at(counts, (ci, li), 1.0)
    tot = counts.sum(1)
    scored = tot > 0
    C = np.full(n_cells, np.nan)
    C[scored] = 1.0 - counts[scored].max(1) / tot[scored]
    return float(np.nanmean(C[scored])) if scored.any() else np.nan, C, int(scored.sum())


def _lineage_composition(cell_index, gene_index, gene_lineage, n_cells):
    """Per-cell lineage-count matrix (n_cells x K) and totals, from the mapping."""
    lut = gene_lineage[gene_index]
    valid = (cell_index >= 0) & (lut >= 0)
    ci, li = cell_index[valid], lut[valid]
    K = int(gene_lineage[gene_lineage >= 0].max()) + 1 if (gene_lineage >= 0).any() else 0
    counts = np.zeros((n_cells, max(K, 1)))
    if K:
        np.add.at(counts, (ci, li), 1.0)
    return counts, counts.sum(1), K


def derive_gene_lineage_from_clusters(expr_cells, cluster_labels, n_genes,
                                     min_cells=5):
    """DATA-DRIVEN, tissue-agnostic lineage map: treat each cell cluster as a
    'lineage' and assign each gene to the cluster in which its mean expression is
    highest. Returns (gene_lineage[n_genes] in [0,K), with -1 for genes not
    confidently enriched). Requires per-CELL expression (cells x genes) and a
    per-cell cluster label. No external panel or tissue knowledge needed."""
    expr = np.asarray(expr_cells, float)
    labels = np.asarray(cluster_labels)
    uniq = [c for c in np.unique(labels) if np.sum(labels == c) >= min_cells]
    if len(uniq) < 2:
        return np.full(n_genes, -1, int)
    # mean expression per cluster (K x G)
    prof = np.vstack([expr[labels == c].mean(0) for c in uniq])
    # z-score each gene across clusters; assign gene to its argmax cluster only if
    # that cluster is a clear peak (z > 0), else leave unassigned (-1).
    mu = prof.mean(0, keepdims=True); sd = prof.std(0, keepdims=True) + 1e-8
    z = (prof - mu) / sd
    best = z.argmax(0); bestz = z.max(0)
    gene_lineage = np.where(bestz > 0.0, best, -1).astype(int)
    return gene_lineage


def validate_refinement(cell_before, cell_after, gene_index, gene_lineage,
                        n_cells, expr_area=None, expr_ecc=None):
    """Multi-criterion before/after validation from the transcript->cell->lineage
    mapping. Returns a dict of paired metrics. Morphology is preserved by
    construction (refinement remaps transcripts only), reported explicitly.
    Metrics requiring downstream re-clustering are out of scope here (Stage 1)."""
    gl = np.asarray(gene_lineage, int)

    def _metrics(ci):
        counts, tot, K = _lineage_composition(ci, gene_index, gl, n_cells)
        scored = tot > 0
        if not scored.any() or K == 0:
            return {"contamination": np.nan, "purity": np.nan,
                    "n_mixed": np.nan, "marker_specificity": np.nan,
                    "n_cells_scored": 0}
        dom = counts[scored].max(1) / tot[scored]
        contamination = float(np.mean(1.0 - dom))
        purity = float(np.mean(dom))
        n_mixed = int(np.sum(dom < 0.5))                     # dominant lineage <50%
        # marker specificity: fraction of each lineage's transcripts that land in
        # cells whose dominant lineage IS that lineage (higher = cleaner markers)
        dom_lineage = np.full(n_cells, -1)
        dom_lineage[scored] = counts[scored].argmax(1)
        lut = gl[gene_index]; v = (ci >= 0) & (lut >= 0)
        cor = 0; tot_l = 0
        for k in range(K):
            mk = v & (lut == k)
            if mk.sum() == 0:
                continue
            in_correct = dom_lineage[ci[mk]] == k
            cor += int(in_correct.sum()); tot_l += int(mk.sum())
        spec = float(cor / tot_l) if tot_l else np.nan
        return {"contamination": contamination, "purity": purity,
                "n_mixed": n_mixed, "marker_specificity": spec,
                "n_cells_scored": int(scored.sum())}

    before = _metrics(cell_before)
    after = _metrics(cell_after)
    moved = int(np.sum((cell_after != cell_before) & (cell_before >= 0)))
    frac_moved = float(moved / max(int(np.sum(cell_before >= 0)), 1))

    out = {
        "contamination_before": before["contamination"],
        "contamination_after": after["contamination"],
        "delta_contamination": (before["contamination"] - after["contamination"]),
        "purity_before": before["purity"],
        "purity_after": after["purity"],
        "delta_purity": (after["purity"] - before["purity"]),
        "n_mixed_before": before["n_mixed"],
        "n_mixed_after": after["n_mixed"],
        "delta_n_mixed": (before["n_mixed"] - after["n_mixed"]),
        "marker_specificity_before": before["marker_specificity"],
        "marker_specificity_after": after["marker_specificity"],
        "delta_marker_specificity": (after["marker_specificity"] - before["marker_specificity"]),
        "n_cells_scored": after["n_cells_scored"],
        "frac_transcripts_reassigned": frac_moved,
        "morphology_preserved": True,   # refinement remaps transcripts only; masks/centroids unchanged
        "note_downstream": ("clustering-stability and neighborhood-consistency are "
                            "out of scope for the side-analysis pass (refined matrix "
                            "does not feed clustering here); evaluate at Stage-2 promotion"),
    }
    dC = out["delta_contamination"]; dS = out["delta_marker_specificity"]
    out["beneficial"] = bool(np.isfinite(dC) and dC > 0 and
                             (not np.isfinite(dS) or dS >= 0))
    return out


# --------------------------------------------------------------------------- #
# core: iterative boundary-transcript reassignment
# --------------------------------------------------------------------------- #
def refine_boundaries(tx_xy: np.ndarray, gene_index: np.ndarray,
                      cell_index0: np.ndarray, centroids_xy: np.ndarray,
                      n_genes: int, tau: float = 0.2, lam: float = 1.0,
                      max_iter: int = 5, tol: float = 1e-3, verbose: bool = True):
    """Refine the transcript->cell assignment at boundaries.

    tx_xy        : (T,2) transcript coordinates (same units as centroids).
    gene_index   : (T,) gene id per transcript in [0, n_genes).
    cell_index0  : (T,) initial cell row per transcript (>=0), or -1 if unassigned.
    centroids_xy : (n_cells,2) nucleus centroids (x,y).
    tau          : ambiguity threshold for the boundary zone.
    lam          : spatial-proximity weight in the score.
    Returns (cell_index_refined, info)."""
    T = tx_xy.shape[0]
    n_cells = centroids_xy.shape[0]
    cell_index = cell_index0.copy()

    i1, d1, i2, d2 = _nearest_two_nuclei(tx_xy, centroids_xy)
    d1s = np.maximum(d1, 1e-8)
    boundary = ((d2 - d1) / d1s) < float(tau)            # ambiguous transcripts
    boundary &= (cell_index0 >= 0)                       # only refine assigned ones
    b_idx = np.where(boundary)[0]
    if verbose:
        print(f"[REFINE] boundary transcripts: {b_idx.size} / {T} "
              f"({100*b_idx.size/max(T,1):.1f}%), tau={tau}, lambda={lam}", flush=True)
    if b_idx.size == 0:
        return cell_index, {"n_boundary": 0, "iterations": 0, "n_reassigned": 0}

    # confidently-assigned = non-boundary transcripts; they anchor the profiles
    conf_mask = (cell_index0 >= 0) & (~boundary)

    def onehot(rows):
        return coo_matrix((np.ones(rows.size), (np.arange(rows.size), gene_index[rows])),
                          shape=(rows.size, n_genes)).tocsr()

    total_reassigned = 0
    it = 0
    for it in range(1, int(max_iter) + 1):
        # (1) profiles from currently confident + already-committed assignments
        conf_rows = np.where(cell_index >= 0)[0]
        theta = _profiles(onehot(conf_rows), cell_index[conf_rows], n_cells)
        logtheta = np.log(theta)

        # (2) reassign each boundary transcript to argmax over its {i1,i2} candidates
        #     (restrict to the two nearest nuclei: local, cheap, interpretable)
        prev = cell_index[b_idx].copy()
        cand = np.stack([i1[b_idx], i2[b_idx]], axis=1)          # (B,2) candidate cells
        g = gene_index[b_idx]
        # transcriptomic term log theta_{cand,g}
        lt = np.stack([logtheta[cand[:, 0], g], logtheta[cand[:, 1], g]], axis=1)
        # spatial term -lambda * distance to each candidate nucleus
        dist = np.stack([d1[b_idx], d2[b_idx]], axis=1)          # d to i1, i2
        score = lt - lam * dist
        pick = score.argmax(1)
        newassign = np.where(pick == 0, cand[:, 0], cand[:, 1])
        cell_index[b_idx] = newassign

        n_changed = int((newassign != prev).sum())
        total_reassigned += n_changed
        if verbose:
            print(f"[REFINE] iter {it}: reassigned {n_changed}", flush=True)
        if n_changed / max(b_idx.size, 1) < tol:
            break

    return cell_index, {"n_boundary": int(b_idx.size), "iterations": it,
                        "n_reassigned": int(total_reassigned)}


# --------------------------------------------------------------------------- #
# pipeline driver
# --------------------------------------------------------------------------- #
def run_refinement(tx_assigned, x_col, y_col, gene_col, cell_col,
                   centroids_xy, genes, out_dir,
                   gene_lineage=None, tau=0.2, lam=1.0, max_iter=5,
                   morphology=None, verbose=True):
    """Refine boundary assignments and write a refined cell-by-gene matrix + a
    before/after validation report. Non-fatal. Never overwrites originals.

    tx_assigned  : DataFrame of assigned transcripts (cols: x, y, gene, cell_id).
    cell_col     : column with the ORIGINAL cell id (1-based label) per transcript.
    centroids_xy : (n_cells,2) centroids aligned to cell rows 0..n_cells-1.
    genes        : list of gene names (len G).
    gene_lineage : (G,) int lineage id per gene (or None to skip contamination)."""
    import pandas as pd
    try:
        os.makedirs(out_dir, exist_ok=True)
        n_cells = centroids_xy.shape[0]
        G = len(genes)

        tx_xy = tx_assigned[[x_col, y_col]].to_numpy(float)
        gmap = {g: k for k, g in enumerate(genes)}
        gene_index = tx_assigned[gene_col].map(gmap).to_numpy()
        ok = ~np.isnan(gene_index)
        gene_index = np.where(ok, gene_index, -1).astype(int)

        # ORIGINAL assignment as 0-based cell rows (assumes cell_col is 1-based label
        # matching centroid row order; -1 for unmapped)
        cell_label = tx_assigned[cell_col].to_numpy()
        cell_index0 = (cell_label - 1).astype(int)
        cell_index0[(cell_index0 < 0) | (cell_index0 >= n_cells) | (~ok)] = -1

        # ---- refine ----
        cell_index1, info = refine_boundaries(
            tx_xy, gene_index, cell_index0, centroids_xy, n_genes=G,
            tau=tau, lam=lam, max_iter=max_iter, verbose=verbose)

        # ---- refined cell-by-gene matrix ----
        def cbg(ci):
            v = (ci >= 0) & (gene_index >= 0)
            return coo_matrix((np.ones(int(v.sum())), (ci[v], gene_index[v])),
                              shape=(n_cells, G)).tocsr()
        M_after = cbg(cell_index1)
        pd.DataFrame(M_after.toarray(), columns=genes).to_csv(
            os.path.join(out_dir, "cell_by_gene_refined.csv"), index_label="cell_row")

        report = {"tau": tau, "lambda": lam, **info,
                  "note": ("transcriptome-guided boundary refinement (Fig. 1a); "
                           "off by default, does not alter masks/centroids; refined "
                           "matrix written alongside the original")}

        # ---- multi-criterion validation: before vs after ----
        if gene_lineage is not None:
            val = validate_refinement(cell_index0, cell_index1, gene_index,
                                      np.asarray(gene_lineage, int), n_cells)
            report.update(val)
            if verbose:
                print(f"[REFINE] deltaC={val['delta_contamination']:+.4f} "
                      f"deltaPurity={val['delta_purity']:+.4f} "
                      f"deltaSpec={val['delta_marker_specificity']:+.4f} "
                      f"mixed {val['n_mixed_before']}->{val['n_mixed_after']} "
                      f"reassigned={val['frac_transcripts_reassigned']:.3f} "
                      f"beneficial={val['beneficial']}", flush=True)

        with open(os.path.join(out_dir, "refinement_report.json"), "w") as f:
            json.dump(report, f, indent=2)
        pd.DataFrame({"transcript_row": np.arange(tx_xy.shape[0]),
                      "cell_before": cell_index0, "cell_after": cell_index1}).to_csv(
            os.path.join(out_dir, "transcript_reassignment.csv"), index=False)
        print(f"[REFINE] wrote refined matrix + report to {out_dir}", flush=True)
        return {"cell_index_refined": cell_index1, "M_after": M_after, "report": report}
    except Exception as e:
        print(f"[REFINE] skipped (non-fatal): {type(e).__name__}: {e}", flush=True)
        return None

