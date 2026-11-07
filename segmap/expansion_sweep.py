#!/usr/bin/env python3
"""
Expansion-radius sweep for SegMap.

Runs the pipeline at several --expand_pixels values (reusing the cached NUCLEUS
mask, so each run is fast and re-segmentation never happens), then collects per-run
quality metrics and writes a summary table + a curve. Purpose: back the choice of
expansion radius with a specificity/quality curve instead of a single number, and
show that ~70 px (10x's ~15 um cell size) is a principled pick rather than
cherry-picked.

Why this is the right experiment: too little expansion under-counts transcripts
(weak power); too much pulls in neighbouring cells' transcripts (specificity drops).
The optimum is where per-cell counts are high but marker specificity / cluster
separability has not yet started to fall.

USAGE:
  1. Edit CONFIG below (XENIUM_DIR and NUCLEUS_MASK especially).
  2. python expansion_sweep.py
  3. Inspect sweep_summary.csv and sweep_curve.png.

Note: this drives the real pipeline via subprocess; it does not re-implement it.
Each run writes to its own --outdir so results don't overwrite each other.
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd

# ============================ CONFIG ============================
XENIUM_DIR = "XeniumData/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs"

# Path to the cached PRE-EXPANSION nucleus mask from a completed run. Reusing this
# explicit path (not 'auto') lets every sweep run share one nucleus segmentation
# while writing to its own outdir.
NUCLEUS_MASK = os.path.join(
    XENIUM_DIR, "cellpose_pipeline_outputs", "data", "full_mask_nuclei.tif")

EXPAND_VALUES = [30, 50, 70, 90, 110]          # pixels; 70 ~= 15 um (10x default)
OUTDIR_BASE = "expansion_sweep"                 # per-run outputs go under here

# Base flags shared by every run (segmentation flags are irrelevant in reuse mode).
BASE_FLAGS = [
    "--gpu", "--no_tune_resolution", "--resolution", "0.3", "--n_pcs", "20",
    "--min_genes_per_cell", "10", "--min_counts_per_cell", "25",
    "--min_cells_per_gene", "3", "--top_gene_cap", "5050",
]
RUN_PIPELINE = True   # set False to only (re)aggregate existing sweep outputs
# ===============================================================


def run_one(expand: int) -> str:
    """Run the pipeline for one expansion value; return its outdir."""
    outdir = os.path.join(OUTDIR_BASE, f"exp{expand}")
    os.makedirs(outdir, exist_ok=True)
    cmd = [sys.executable, "-m", "segmap.run_segmap",
           "--xenium_dir", XENIUM_DIR,
           "--reuse_mask", NUCLEUS_MASK,
           "--outdir", outdir,
           "--expand_pixels", str(expand)] + BASE_FLAGS
    print("\n>>> RUN expand=%d\n    %s" % (expand, " ".join(cmd)), flush=True)
    subprocess.run(cmd, check=True)
    return outdir


def _data_dir(outdir: str) -> str:
    # The pipeline nests outputs under <outdir>/data (or <outdir> directly).
    for cand in (os.path.join(outdir, "data"), outdir):
        if os.path.isdir(cand):
            return cand
    return outdir


def collect_metrics(outdir: str, expand: int) -> dict:
    """Pull quality metrics from one run's outputs. Robust to missing files."""
    d = _data_dir(outdir)
    row = {"expand_px": expand, "expand_um": round(expand * 0.2125, 2)}

    # --- separability + cell/cluster counts ---
    sep_path = os.path.join(d, "cluster_separability.csv")
    try:
        sep = pd.read_csv(sep_path)
        n = sep["n"].to_numpy(dtype=float)
        sil = sep["mean_silhouette"].to_numpy(dtype=float)
        row["n_cells"] = int(n.sum())
        row["n_clusters"] = int(len(sep))
        row["weighted_silhouette"] = float(np.average(sil, weights=n)) if n.sum() else np.nan
    except Exception as e:
        print(f"  [warn] separability read failed ({e})")

    # --- median transcripts per cell ---
    for name in ("cell_by_gene_for_umap_final.csv", "cell_by_gene.csv"):
        cg_path = os.path.join(d, name)
        if os.path.exists(cg_path):
            try:
                cg = pd.read_csv(cg_path)
                counts = cg.iloc[:, 1:].to_numpy(dtype=float).sum(axis=1)
                row["median_counts_per_cell"] = float(np.median(counts))
                row["mean_counts_per_cell"] = float(np.mean(counts))
            except Exception as e:
                print(f"  [warn] cell_by_gene read failed ({e})")
            break

    # --- marker strength + specificity (if detection columns exist) ---
    mk_path = os.path.join(d, "cluster_markers_all_cells_final.csv")
    try:
        mk = pd.read_csv(mk_path)
        fdr = "pvals_adj" if "pvals_adj" in mk.columns else None
        lfc = "logfoldchanges" if "logfoldchanges" in mk.columns else None
        sig = mk
        if fdr and lfc:
            sig = mk[(mk[fdr] < 0.05) & (mk[lfc] > 0)]
        if lfc:
            row["median_log2FC"] = float(np.median(sig[lfc])) if len(sig) else np.nan
        # detection-fraction specificity if present
        p1 = next((c for c in ("pct_nz_group", "pct.1", "pct1") if c in mk.columns), None)
        p2 = next((c for c in ("pct_nz_reference", "pct.2", "pct2") if c in mk.columns), None)
        if p1 and p2:
            row["mean_specificity"] = float(np.mean(sig[p1] - sig[p2])) if len(sig) else np.nan
        else:
            row["mean_specificity"] = np.nan  # not exported; see note in the summary
    except Exception as e:
        print(f"  [warn] marker read failed ({e})")

    return row


def make_plot(df: pd.DataFrame, path: str = "sweep_curve.png") -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"[plot] matplotlib unavailable ({e}); skipping plot.")
        return
    x = df["expand_px"]
    fig, ax1 = plt.subplots(figsize=(7, 5))
    ax1.set_xlabel("expansion radius (px)   [~%.2f um/px]" % 0.2125)
    l1 = ax1.plot(x, df.get("median_counts_per_cell"), "o-", color="tab:blue",
                  label="median counts/cell")
    ax1.set_ylabel("median transcripts / cell", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax1.axvline(70, ls="--", color="gray", alpha=0.6)
    ax1.annotate("~15 um (10x)", (70, ax1.get_ylim()[1]), fontsize=8, color="gray")

    ax2 = ax1.twinx()
    ycol = "mean_specificity" if df["mean_specificity"].notna().any() else "weighted_silhouette"
    l2 = ax2.plot(x, df.get(ycol), "s-", color="tab:red", label=ycol)
    ax2.set_ylabel(ycol, color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")

    lines = l1 + l2
    ax1.legend(lines, [ln.get_label() for ln in lines], loc="best", fontsize=8)
    plt.title("Expansion sweep: counts vs quality (pick where quality peaks)")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    print(f"[plot] wrote {path}")


def main():
    rows = []
    for e in EXPAND_VALUES:
        outdir = os.path.join(OUTDIR_BASE, f"exp{e}")
        if RUN_PIPELINE:
            try:
                outdir = run_one(e)
            except subprocess.CalledProcessError as err:
                print(f"[error] pipeline failed for expand={e}: {err}")
                continue
        rows.append(collect_metrics(outdir, e))

    df = pd.DataFrame(rows).sort_values("expand_px")
    df.to_csv("sweep_summary.csv", index=False)
    print("\n===== EXPANSION SWEEP SUMMARY =====")
    print(df.to_string(index=False))
    make_plot(df)
    print("\nPick the radius where transcripts/cell is high but "
          "specificity/silhouette has NOT started to drop — that is the point past "
          "which expansion begins capturing neighbours' transcripts.")
    if "mean_specificity" in df and df["mean_specificity"].isna().all():
        print("\nNOTE: marker specificity (pct.1-pct.2) was not in the marker CSV, so "
              "the quality curve uses cluster silhouette instead. To get true "
              "specificity, export detection fractions in the marker table.")


if __name__ == "__main__":
    main()
