#!/usr/bin/env python3
"""
segmentation_sensitivity.py -- sensitivity analysis for the SegMap
segmentation / cell-reconstruction / spatial-network hyperparameters.

DEFAULT MODE = geometry-only sweep (mask + transcripts; no downstream clustering):
  reconstructs cells from the saved NUCLEUS MASK under each hyperparameter setting
  and reports, per setting per dataset:
      n_cells, median_cell_area_px, median_solidity, frac_transcripts_assigned
  and the Jaccard overlap of the cell FOOTPRINT vs the default (how much the set
  of assigned pixels changes) as a mask-level concordance metric.

OPTIONAL MODE = clustering ARI (needs YOUR assignment+clustering functions):
  set USE_CLUSTERING=True and complete cluster_from_labels() to additionally
  report the adjusted Rand index of the downstream clustering vs the default.

Hyperparameters swept (from the Methods):
  r            : nucleus-expansion radius (px)     -- Tier 1, most consequential
  min_cell_area / max_cell_area / min_solidity     -- geometry QC filters
  (k, the spatial-graph neighbours, only affects downstream clustering, so it is
   swept only in the optional clustering mode.)

Calls SegMap's own resolve_from_nucleus_mask(); does not reimplement it.

Deps: numpy, pandas, scikit-image, matplotlib, tifffile, pyarrow
      (+ scikit-learn and the segmap package only for the optional ARI mode).
"""
from __future__ import annotations
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# EDIT THIS BLOCK
# --------------------------------------------------------------------------- #
PIXEL_SIZE_UM = 0.2125          # Xenium pixel size (um/px); um -> px = value / PIXEL_SIZE_UM

# Raw Xenium transcripts per dataset NAME (from your `find` output).
TRANSCRIPTS = {
    "AD":      "XeniumData/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/transcripts.parquet",
    "GBM":     "XeniumData/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon_outs/transcripts.parquet",
    "HEALTHY": "XeniumData/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon_outs/transcripts.parquet",
    "LUNG":    "XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs/transcripts.parquet",
}

# Candidate column names in transcripts.parquet (first present is used).
X_COLS = ["x_location", "x", "x_centroid", "x_pixel"]
Y_COLS = ["y_location", "y", "y_centroid", "y_pixel"]
# True if the x/y columns are in MICRONS (Xenium default) -> divide by PIXEL_SIZE_UM.
TRANSCRIPT_COORDS_IN_MICRONS = True

DEFAULTS = dict(
    r=15,                 # <<< SET to the expansion radius used in the paper
    min_cell_area=15, max_cell_area=0, min_solidity=0.0,
    resolve_mode="expand_labels",
    k=15,                 # only used in the optional clustering mode
)
GRIDS = {
    "r":             [8, 10, 12, 15, 18, 22, 26, 30],
    "min_cell_area": [5, 15, 30, 50],
    "max_cell_area": [0, 2000, 4000],
    "min_solidity":  [0.0, 0.8, 0.9],
}
PLABEL = {"r": "expansion radius r (px)", "min_cell_area": "min cell area (px)",
          "max_cell_area": "max cell area (px)", "min_solidity": "min solidity",
          "k": "kNN neighbours k"}

USE_CLUSTERING = False          # set True after wiring cluster_from_labels()


# --------------------------------------------------------------------------- #
# IO
# --------------------------------------------------------------------------- #
def _pick(cols, candidates, what):
    for c in candidates:
        if c in cols:
            return c
    raise KeyError(f"none of {candidates} found for {what}; columns are {list(cols)}")

def load_segmap_state(path, name):
    import tifffile
    mask = tifffile.imread(os.path.join(path, "full_mask_nuclei.tif")).astype(np.int32)
    tx = pd.read_parquet(TRANSCRIPTS[name])
    xc = _pick(tx.columns, X_COLS, "x"); yc = _pick(tx.columns, Y_COLS, "y")
    x = tx[xc].to_numpy(float); y = tx[yc].to_numpy(float)
    if TRANSCRIPT_COORDS_IN_MICRONS:
        x = x / PIXEL_SIZE_UM; y = y / PIXEL_SIZE_UM
    print(f"    [{name}] mask {mask.shape}, {len(tx):,} transcripts "
          f"(x='{xc}', y='{yc}', microns={TRANSCRIPT_COORDS_IN_MICRONS})", flush=True)
    return {"nucleus_mask": mask, "tx_px": np.rint(np.c_[x, y]).astype(int),
            "shape": mask.shape, "name": name}


# --------------------------------------------------------------------------- #
# One reconstruction under a given hyperparameter setting
# --------------------------------------------------------------------------- #
def run_from_mask(seg, params):
    from segmap.resolve_overlaps import resolve_from_nucleus_mask
    res = resolve_from_nucleus_mask(
        seg["nucleus_mask"], expand_pixels=int(params["r"]),
        mode=params["resolve_mode"],
        min_area=int(params["min_cell_area"]),
        max_area=(int(params["max_cell_area"]) or None),
        min_solidity=float(params["min_solidity"]), verbose=False)
    labels_img = res["labels"]; morph = res["morphology"] or {}

    # transcript assignment fraction: transcripts landing on a labelled pixel
    H, W = labels_img.shape
    xy = seg["tx_px"]; xs, ys = xy[:, 0], xy[:, 1]
    inb = (xs >= 0) & (xs < W) & (ys >= 0) & (ys < H)
    lab = np.zeros(xs.shape[0], dtype=labels_img.dtype)
    lab[inb] = labels_img[ys[inb], xs[inb]]
    frac = float((lab > 0).mean())

    areas = np.asarray(morph.get("area", []), float)
    solid = np.asarray(morph.get("solidity", []), float)
    out = dict(
        n_cells=int(labels_img.max()),
        median_cell_area_px=float(np.median(areas)) if areas.size else np.nan,
        median_solidity=float(np.median(solid)) if solid.size else np.nan,
        frac_transcripts_assigned=frac,
        foreground=(labels_img > 0),          # for footprint Jaccard vs default
    )
    if USE_CLUSTERING:
        out.update(cluster_from_labels(seg, labels_img, params))   # see stub below
    return out


def footprint_jaccard(fg_a, fg_b):
    inter = np.logical_and(fg_a, fg_b).sum()
    union = np.logical_or(fg_a, fg_b).sum()
    return 1.0 if union == 0 else float(inter / union)


# --------------------------------------------------------------------------- #
# Sweep one dataset
# --------------------------------------------------------------------------- #
def run_dataset(name, seg, grids):
    print(f"[{name}] default reconstruction (r={DEFAULTS['r']}) ...", flush=True)
    default = run_from_mask(seg, DEFAULTS)
    rows = []
    for param, values in grids.items():
        if param == "k" and not USE_CLUSTERING:
            continue
        for v in values:
            params = dict(DEFAULTS); params[param] = v
            pert = default if params == DEFAULTS else run_from_mask(seg, params)
            row = dict(dataset=name, parameter=param, label=PLABEL[param], value=v,
                       is_default=(params == DEFAULTS),
                       n_cells=pert["n_cells"],
                       median_cell_area_px=pert["median_cell_area_px"],
                       median_solidity=pert["median_solidity"],
                       frac_transcripts_assigned=pert["frac_transcripts_assigned"],
                       footprint_jaccard=footprint_jaccard(default["foreground"],
                                                           pert["foreground"]))
            if USE_CLUSTERING and "labels" in pert and "labels" in default:
                from sklearn.metrics import adjusted_rand_score
                ids = np.intersect1d(default["cell_ids"], pert["cell_ids"])
                if ids.size >= 2:
                    dm = dict(zip(default["cell_ids"], default["labels"]))
                    pm = dict(zip(pert["cell_ids"], pert["labels"]))
                    row["ARI_clustering"] = float(adjusted_rand_score(
                        [dm[i] for i in ids], [pm[i] for i in ids]))
            rows.append(row)
            print(f"    [{name}] {param}={v}: n_cells={row['n_cells']} "
                  f"area={row['median_cell_area_px']:.0f} "
                  f"fracTx={row['frac_transcripts_assigned']:.3f} "
                  f"footJ={row['footprint_jaccard']:.3f}", flush=True)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def make_figure(df, params, outpath):
    metrics = [("n_cells", "cell count"),
               ("median_cell_area_px", "median cell area (px)"),
               ("frac_transcripts_assigned", "fraction transcripts assigned"),
               ("footprint_jaccard", "footprint Jaccard vs default")]
    if "ARI_clustering" in df.columns:
        metrics.insert(0, ("ARI_clustering", "ARI vs default (clustering)"))
    metrics = [m for m in metrics if m[0] in df.columns]
    params = [p for p in params if p in df.parameter.unique()]
    fig, axes = plt.subplots(len(metrics), len(params),
                             figsize=(3.0 * len(params), 2.2 * len(metrics)), squeeze=False)
    ds = sorted(df.dataset.unique()); cmap = plt.get_cmap("tab10")
    for r_, (mk, ml) in enumerate(metrics):
        for c_, p in enumerate(params):
            ax = axes[r_][c_]; sub = df[df.parameter == p]
            for di, dn in enumerate(ds):
                g = sub[sub.dataset == dn]
                try: x = g.value.astype(float).to_numpy()
                except (TypeError, ValueError): x = np.arange(len(g))
                o = np.argsort(x)
                ax.plot(np.asarray(x)[o], g[mk].to_numpy()[o], "-o", ms=4,
                        color=cmap(di % 10), label=dn)
            if r_ == 0: ax.set_title(PLABEL[p], fontsize=8)
            if c_ == 0: ax.set_ylabel(ml, fontsize=8)
            if mk in ("footprint_jaccard", "ARI_clustering", "frac_transcripts_assigned"):
                ax.set_ylim(0, 1.02)
            ax.tick_params(labelsize=7)
    axes[0][-1].legend(fontsize=7, frameon=False)
    fig.suptitle("Sensitivity of SegMap cell reconstruction to its hyperparameters", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, bbox_inches="tight")
    print(f"[fig] wrote {outpath}", flush=True)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", required=True, help="name=path pairs")
    ap.add_argument("--outdir", default="seg_sensitivity_out")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    all_df = []
    for spec in args.datasets:
        name, path = spec.split("=", 1)
        seg = load_segmap_state(path, name)
        all_df.append(run_dataset(name, seg, GRIDS))
    df = pd.concat(all_df, ignore_index=True)

    df.to_csv(os.path.join(args.outdir, "segmentation_sensitivity.csv"), index=False)
    keep = ["dataset", "parameter", "label", "value", "is_default", "n_cells",
            "median_cell_area_px", "median_solidity", "frac_transcripts_assigned",
            "footprint_jaccard"] + (["ARI_clustering"] if "ARI_clustering" in df.columns else [])
    df[keep].sort_values(["dataset", "parameter", "value"]).to_csv(
        os.path.join(args.outdir, "SupplementaryTable_segmentation_sensitivity.csv"), index=False)
    make_figure(df, list(GRIDS.keys()),
                os.path.join(args.outdir, "SupplementaryFigure_segmentation_sensitivity.pdf"))
    print(f"[done] wrote table + figure to {args.outdir}/", flush=True)


# --------------------------------------------------------------------------- #
# OPTIONAL clustering hook -- complete to enable ARI (set USE_CLUSTERING=True).
# Must return {"cell_ids": <array>, "labels": <array>} for the reconstructed cells.
# --------------------------------------------------------------------------- #
def cluster_from_labels(seg, labels_img, params):
    """Assign transcripts to `labels_img`, build the k-NN graph (params['k']),
    cluster, and return per-cell labels. Wire to YOUR pipeline functions, e.g.:

        from segmap.pipeline import assign_transcripts_to_cells_fullres
        from segmap.preprocessing import exact_kmeans_clustering
        adata = assign_transcripts_to_cells_fullres(labels_img, seg["transcripts"])
        adata = exact_kmeans_clustering(adata, n_clusters=..., n_pcs=..., seed=0)
        return {"cell_ids": np.asarray(adata.obs_names),
                "labels":   adata.obs["cluster"].to_numpy()}
    """
    raise NotImplementedError("Complete cluster_from_labels() to enable ARI mode.")


if __name__ == "__main__":
    main()
