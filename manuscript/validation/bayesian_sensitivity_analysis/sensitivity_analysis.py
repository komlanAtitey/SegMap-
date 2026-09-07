#!/usr/bin/env python3
"""
sensitivity_analysis.py -- hyperparameter sensitivity analysis + visualization
for the SegMap graph-spectral Bayesian marker model (bayes_markers.py).

For each tunable setting (k, m, d, alpha, Laplacian normalization) the model is
refit around its default, and the perturbed fit is compared with the default by
three concordance metrics:
    (i)   Spearman rank correlation of per-gene posteriors p_j = p(X_j=1 | y_j);
    (ii)  Jaccard overlap of the Bayesian-FDR marker sets (nominal FDR < alpha);
    (iii) change in the number of selected markers.

Outputs (under --outdir):
    sensitivity_concordance.csv            tidy table, one row per setting/dataset
    SupplementaryTable_sensitivity.csv     same, formatted for the supplement
    SupplementaryFigure_sensitivity.pdf    concordance curves per parameter

This calls the ACTUAL model functions in bayes_markers.py; it does not
reimplement them. Point --datasets at your normalized expression + centroids.

Deps: numpy, pandas, scipy, matplotlib (+ the segmap package on PYTHONPATH).
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- import the real model -------------------------------------------------- #
# Adjust if your package layout differs (e.g. `from segmap.bayes_markers import ...`).
try:
    from bayes_markers import build_spatial_basis, fit_all_genes, select_markers_bfdr
except Exception:
    from segmap.bayes_markers import build_spatial_basis, fit_all_genes, select_markers_bfdr


# =========================================================================== #
# Defaults and perturbation grids  (match Supplementary Information X)
# =========================================================================== #
DEFAULTS = dict(k=15, m=100, d=1, alpha=0.05, normalized=True)

GRIDS = {
    "k":          [5, 10, 15, 20, 30],
    "m":          [50, 100, 150, 200],
    "d":          [1, 2],
    "alpha":      [0.01, 0.05, 0.10],
    "normalized": [True, False],
}
PARAM_LABEL = {
    "k": "kNN neighbours k", "m": "spectral modes m",
    "d": "BIC penalty d", "alpha": "FDR level \u03b1",
    "normalized": "Laplacian (norm. vs comb.)",
}


# =========================================================================== #
# One model fit for a given hyperparameter setting
#   -> DataFrame(gene, p_j) and the selected-marker gene set
# =========================================================================== #
def fit_once(expr, coords, genes, k, m, d, alpha, normalized):
    """Run the real model once and return (per-gene posteriors, marker set)."""
    basis = build_spatial_basis(coords, k=int(k), m=int(m), normalized=bool(normalized))
    fit = fit_all_genes(expr, basis, dpar=int(d))
    p = np.asarray(fit["posterior"], float)
    sel = select_markers_bfdr(p, alpha=float(alpha))          # index array
    markers = set(np.asarray(genes)[sel].tolist())
    return p, markers


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def concordance(p_def, m_def, p_pert, m_pert):
    rho = spearmanr(p_def, p_pert).correlation
    return dict(
        spearman_rho=float(rho),
        jaccard_marker=float(jaccard(m_def, m_pert)),
        n_markers_def=int(len(m_def)),
        n_markers_pert=int(len(m_pert)),
        delta_markers=int(len(m_pert) - len(m_def)),
    )


# =========================================================================== #
# Sweep one dataset
# =========================================================================== #
def run_dataset(name, expr, coords, genes):
    print(f"[{name}] default fit (k={DEFAULTS['k']}, m={DEFAULTS['m']}, "
          f"d={DEFAULTS['d']}, alpha={DEFAULTS['alpha']}) ...", flush=True)
    p_def, m_def = fit_once(expr, coords, genes, **DEFAULTS)

    rows = []
    for param, values in GRIDS.items():
        for v in values:
            setting = dict(DEFAULTS); setting[param] = v
            # skip the exact default (already computed) but keep it as a reference row
            p_pert, m_pert = (p_def, m_def) if setting == DEFAULTS else \
                fit_once(expr, coords, genes, **setting)
            cc = concordance(p_def, m_def, p_pert, m_pert)
            rows.append(dict(dataset=name, parameter=param, label=PARAM_LABEL[param],
                             value=v, is_default=(setting == DEFAULTS), **cc))
            print(f"    [{name}] {param}={v}: rho={cc['spearman_rho']:.3f} "
                  f"J={cc['jaccard_marker']:.3f} dN={cc['delta_markers']:+d}", flush=True)
    return pd.DataFrame(rows)


# =========================================================================== #
# Visualization: concordance curves per parameter
# =========================================================================== #
def make_figure(df, outpath):
    params = list(GRIDS.keys())
    metrics = [("spearman_rho", "Spearman \u03c1 (posteriors)"),
               ("jaccard_marker", "Jaccard (marker sets)")]
    fig, axes = plt.subplots(len(metrics), len(params),
                             figsize=(3.0 * len(params), 5.2), squeeze=False)
    datasets = sorted(df["dataset"].unique())
    cmap = plt.get_cmap("tab10")
    for r, (mkey, mlab) in enumerate(metrics):
        for c, param in enumerate(params):
            ax = axes[r][c]
            sub = df[df["parameter"] == param]
            for di, ds in enumerate(datasets):
                d = sub[sub["dataset"] == ds].copy()
                # numeric x where possible; else categorical positions
                try:
                    x = d["value"].astype(float).to_numpy()
                except (TypeError, ValueError):
                    x = np.arange(len(d))
                order = np.argsort(x)
                ax.plot(np.asarray(x)[order], d[mkey].to_numpy()[order],
                        "-o", ms=4, color=cmap(di % 10), label=ds)
            ax.set_ylim(0, 1.02)
            if r == 0:
                ax.set_title(PARAM_LABEL[param], fontsize=9)
            if c == 0:
                ax.set_ylabel(mlab, fontsize=9)
            ax.axhline(0.9, ls=":", lw=0.8, color="grey")
            ax.tick_params(labelsize=8)
    axes[0][-1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("Sensitivity of the Bayesian spatial marker model to its hyperparameters",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, bbox_inches="tight")
    print(f"[fig] wrote {outpath}", flush=True)


# =========================================================================== #
# Dataset loading  --  EDIT to your file layout
# =========================================================================== #
def load_dataset(path):
    """Read SegMap outputs: cell-by-gene matrix + centroid coordinates.
    `path` is the dataset's data dir (…/cellpose_pipeline_outputs_X/data)."""
    import numpy as np, pandas as pd, os

    cg = pd.read_csv(os.path.join(path, "cell_by_gene_for_umap_final.csv"), index_col=0)
    sp = pd.read_csv(os.path.join(path, "spatial_clusters_final.csv"))

    # align cells present in both tables (coords keyed by the same cell id)
    # -- adjust the id column name if yours differs:
    id_col = sp.columns[0] if "cell" not in sp.columns[0].lower() else sp.columns[0]
    sp = sp.set_index(sp.columns[0])
    common = cg.index.intersection(sp.index)
    cg, sp = cg.loc[common], sp.loc[common]

    genes  = list(cg.columns)
    counts = cg.to_numpy(float)                       # raw counts, cells x genes
    coords = sp[["x_centroid_px", "y_centroid_px"]].to_numpy(float)

    # normalize exactly as the clustering pipeline does: CP10k + log1p
    lib = counts.sum(1, keepdims=True); lib[lib == 0] = 1.0
    expr = np.log1p(counts / lib * 1e4)
    return expr, coords, genes


# =========================================================================== #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", required=True,
                    help="name=path pairs, e.g. AD=data/AD LUNG=data/LUNG")
    ap.add_argument("--outdir", default="sensitivity_out")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    all_df = []
    for spec in args.datasets:
        name, path = spec.split("=", 1)
        expr, coords, genes = load_dataset(path)
        all_df.append(run_dataset(name, expr, coords, genes))
    df = pd.concat(all_df, ignore_index=True)

    csv1 = os.path.join(args.outdir, "sensitivity_concordance.csv")
    df.to_csv(csv1, index=False)
    supp = df[["dataset", "parameter", "label", "value", "is_default",
               "spearman_rho", "jaccard_marker",
               "n_markers_def", "n_markers_pert", "delta_markers"]] \
        .sort_values(["dataset", "parameter", "value"])
    supp.to_csv(os.path.join(args.outdir, "SupplementaryTable_sensitivity.csv"), index=False)
    make_figure(df, os.path.join(args.outdir, "SupplementaryFigure_sensitivity.pdf"))
    print(f"[done] wrote table + figure to {args.outdir}/", flush=True)


if __name__ == "__main__":
    main()
