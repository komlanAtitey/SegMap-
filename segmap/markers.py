"""Cluster marker detection (fast + Scanpy Wilcoxon) and auto-annotation."""
from .common import *  # noqa: F401,F403


def fast_cluster_markers_foldchange(
    adata: ad.AnnData,
    cluster_key: str = "cluster",
    top_n: int = 10,
    min_cells_in: int = 10,
    min_pct_in: float = 0.10,
    min_pct_diff: float = 0.0,
) -> pd.DataFrame:  # fast, robust markers
    """Fast marker detection with proper normalization + detection-fraction filtering.

    Accuracy upgrades vs. the old version:
      * Counts are CP10k-normalized + log1p BEFORE comparing means, so a cell's
        sequencing depth no longer dominates the fold change.
      * Fold change uses a data-scale pseudocount (log-space difference) instead of
        a 1e-6 pseudocount that made near-zero out-group genes explode to huge FC.
      * Genes must be detected in at least `min_pct_in` of in-cluster cells
        (and optionally exceed the out-group detection rate by `min_pct_diff`),
        which removes the rare-noise genes that previously won the ranking.
      * Ranking is a specificity score (log2FC weighted by detection fraction),
        not raw log2FC, so the top genes are genuinely cluster-defining.
    """
    Xc = adata.X  # raw counts matrix
    Xc = Xc.toarray() if sp.issparse(Xc) else np.asarray(Xc)  # dense
    Xc = Xc.astype(np.float32, copy=False)  # float

    # ---- CP10k normalization + log1p (depth-corrected expression) ----
    lib = Xc.sum(axis=1, keepdims=True)  # per-cell library size
    lib[lib == 0] = 1.0  # avoid divide-by-zero
    Xln = np.log1p((Xc / lib) * 1e4)  # log-normalized expression
    detected = (Xc > 0)  # boolean detection matrix

    genes = np.array(adata.var_names).astype(str)  # gene names
    clusters = adata.obs[cluster_key].astype(str).to_numpy()  # labels
    uniq = np.unique(clusters)  # unique clusters
    out_rows: List[Dict[str, Any]] = []  # output rows

    for cl in uniq:  # loop clusters
        idx_in = (clusters == cl)  # in-cluster mask
        idx_out = ~idx_in  # out-cluster mask
        n_in = int(idx_in.sum())  # in-cluster size
        if n_in < int(min_cells_in):  # too small to trust
            out_rows.append({"cluster": cl, "gene": "(cluster too small)", "log2FC": np.nan,
                             "mean_in": np.nan, "mean_out": np.nan, "pct_in": np.nan,
                             "pct_out": np.nan, "score": np.nan, "n_cells_in": n_in})  # placeholder
            continue  # next cluster

        mean_in = Xln[idx_in].mean(axis=0)  # mean log-expr in
        mean_out = Xln[idx_out].mean(axis=0) if int(idx_out.sum()) > 0 else np.zeros_like(mean_in)  # mean log-expr out
        pct_in = detected[idx_in].mean(axis=0)  # detection fraction in
        pct_out = detected[idx_out].mean(axis=0) if int(idx_out.sum()) > 0 else np.zeros_like(pct_in)  # detection fraction out

        # log2 fold change of expm1(mean) with pseudocount 1 (stable, depth-corrected)
        log2fc = np.log2((np.expm1(mean_in) + 1.0) / (np.expm1(mean_out) + 1.0))  # log2FC

        # specificity score: up-regulation x how cluster-specific the detection is
        score = log2fc * pct_in * (pct_in - pct_out + 1.0)  # combined ranking score

        # detection-fraction gate removes rare-noise genes
        keep = (pct_in >= float(min_pct_in)) & ((pct_in - pct_out) >= float(min_pct_diff)) & (log2fc > 0)  # mask
        if not np.any(keep):  # fall back to FC if nothing passes the gate
            keep = (log2fc > 0)  # relaxed
        if not np.any(keep):  # truly nothing up-regulated
            keep = np.ones_like(log2fc, dtype=bool)  # keep all as last resort

        cand = np.where(keep)[0]  # candidate gene indices
        order = cand[np.argsort(score[cand])[::-1]][: int(top_n)]  # top by score
        for j in order:  # emit rows
            out_rows.append({"cluster": cl, "gene": genes[j], "log2FC": float(log2fc[j]),
                             "mean_in": float(mean_in[j]), "mean_out": float(mean_out[j]),
                             "pct_in": float(pct_in[j]), "pct_out": float(pct_out[j]),
                             "score": float(score[j]), "n_cells_in": n_in})  # row
    return pd.DataFrame(out_rows)  # return dataframe


def auto_annotate_clusters(adata, cluster_key="cluster", top_n=5):
    """Unsupervised cluster annotation using the most cluster-specific genes.

    Prefers Scanpy's Wilcoxon rank_genes_groups (statistically grounded) when
    Scanpy is available and the matrix is large enough; otherwise falls back to
    the robust fast fold-change markers. Placeholder/empty genes are skipped so
    labels are always real gene names.
    """
    cluster_names = {}

    # ---- Preferred path: Scanpy Wilcoxon DE (more accurate labels) ----
    if SCANPY_AVAILABLE and int(adata.n_obs) >= 50:
        try:
            am = adata.copy()
            sc.pp.normalize_total(am, target_sum=1e4)
            sc.pp.log1p(am)
            am.obs[cluster_key] = am.obs[cluster_key].astype("category")
            counts = am.obs[cluster_key].value_counts()
            valid = counts[counts >= 5].index
            am = am[am.obs[cluster_key].isin(valid)].copy()
            if am.obs[cluster_key].nunique() >= 2:
                sc.tl.rank_genes_groups(am, groupby=cluster_key, method="wilcoxon")
                names = am.uns["rank_genes_groups"]["names"]
                for cl in names.dtype.names:
                    genes = [g for g in list(names[cl])[: top_n + 3] if g and not str(g).startswith("(")]
                    cluster_names[str(cl)] = "_".join(genes[:3]) if genes else str(cl)
        except Exception as e:
            print(f"[ANNOTATE] Scanpy annotation failed ({e}); using fast markers.", flush=True)

    # ---- Fallback / fill-in: robust fast markers ----
    if not cluster_names:
        markers = fast_cluster_markers_foldchange(adata, cluster_key=cluster_key, top_n=top_n)
        for cl in markers["cluster"].unique():
            genes = [g for g in markers[markers["cluster"] == cl]["gene"].tolist()
                     if g and not str(g).startswith("(")]
            cluster_names[str(cl)] = "_".join(genes[:3]) if genes else str(cl)

    # ---- Ensure every cluster has a label ----
    for cl in adata.obs[cluster_key].astype(str).unique():
        cluster_names.setdefault(str(cl), str(cl))

    return cluster_names


def compute_scanpy_markers_all_cells(
    adata,
    groupby,
    output_path
):
    """
    Compute cluster markers using all cells (Scanpy Wilcoxon test).

    Outputs:
        output_path
        output_path.replace(".csv", "_top20.csv")
    """

    import scanpy as sc
    import pandas as pd
    import numpy as np

    print(
        f"[MARKERS] Computing markers for '{groupby}'...",
        flush=True
    )

    adata_m = adata.copy()

    # ==========================================================
    # Normalize and log-transform
    # ==========================================================
    sc.pp.normalize_total(
        adata_m,
        target_sum=1e4
    )

    sc.pp.log1p(adata_m)

    # ==========================================================
    # Remove noise cluster (if present)
    # ==========================================================
    if "noise" in adata_m.obs[groupby].astype(str).unique():

        adata_m = adata_m[
            adata_m.obs[groupby].astype(str) != "noise"
        ].copy()

        print(
            "[MARKERS] Removed noise cluster",
            flush=True
        )

    # ==========================================================
    # Remove tiny groups before marker detection
    # ==========================================================
    group_counts = adata_m.obs[groupby].value_counts()

    valid_groups = group_counts[group_counts >= 5].index

    adata_m = adata_m[
        adata_m.obs[groupby].isin(valid_groups)
    ].copy()

    print(
        f"[MARKERS] Retained {len(valid_groups)} groups "
        f"with >=5 cells",
        flush=True
    )

    # Must have at least 2 groups for DE analysis
    if len(valid_groups) < 2:

        print(
            "[MARKERS] Not enough valid groups for DE analysis.",
            flush=True
        )

        return

    # ==========================================================
    # Compute markers using ALL genes
    # ==========================================================
    try:

        sc.tl.rank_genes_groups(
            adata_m,
            groupby=groupby,
            method="wilcoxon"
        )

    except Exception as e:

        print(
            f"[MARKERS] rank_genes_groups failed: {e}",
            flush=True
        )

        return

    # ==========================================================
    # Extract results
    # ==========================================================
    df = sc.get.rank_genes_groups_df(
        adata_m,
        group=None
    )
    
    # ==========================================================
    # Keep standard columns
    # ==========================================================
    keep_cols = [
        "group",
        "names",
        "scores",
        "logfoldchanges",
        "pvals",
        "pvals_adj"
    ]

    df = df[keep_cols].copy()

    # ==========================================================
    # Save complete marker table
    # ==========================================================
    df.to_csv(
        output_path,
        index=False
    )

    print(
        f"[MARKERS] Saved full table: {output_path}",
        flush=True
    )

    # ==========================================================
    # Significant markers only
    # ==========================================================
    sig_df = df[
        (df["pvals_adj"] < 0.05) &
        (df["logfoldchanges"] > 0.5)
    ].copy()

    sig_path = output_path.replace(
        ".csv",
        "_significant.csv"
    )

    sig_df.to_csv(
        sig_path,
        index=False
    )

    print(
        f"[MARKERS] Saved significant markers: {sig_path}",
        flush=True
    )

    # ==========================================================
    # Top 20 markers per cluster
    # ==========================================================
    top20 = (
        sig_df
        .sort_values(
            ["group", "scores"],
            ascending=[True, False]
        )
        .groupby("group")
        .head(20)
        .reset_index(drop=True)
    )

    top20_path = output_path.replace(
        ".csv",
        "_top20.csv"
    )

    top20.to_csv(
        top20_path,
        index=False
    )

    print(
        f"[MARKERS] Saved top20 markers: {top20_path}",
        flush=True
    )

    # ==========================================================
    # Summary
    # ==========================================================
    print(
        "[MARKERS] Marker discovery complete.",
        flush=True
    )
