"""AnnData QC filtering and PCA/HVG embedding."""
from .common import *  # noqa: F401,F403


def safe_filter_adata(adata: ad.AnnData, min_cells_floor: int = 50,
                      min_genes: int = None, min_cells: int = None,
                      min_counts: int = 0) -> ad.AnnData:  # safe filtering
    if int(adata.n_obs) == 0:  # empty
        return adata  # return

    # Total-counts-per-cell floor: removes the very sparse nuclei that don't sit
    # cleanly in any cluster and smear the UMAP, without the all-or-nothing jump
    # that raising min_genes alone forces.
    if min_counts and int(min_counts) > 0:
        tot = np.asarray(adata.X.sum(axis=1)).ravel()
        keep = tot >= int(min_counts)
        n_drop = int((~keep).sum())
        a_c = adata[keep].copy()
        if int(a_c.n_obs) > 0:
            print(f"[FILTER] min_counts_per_cell={int(min_counts)} dropped {n_drop} cells "
                  f"-> {a_c.n_obs} remain", flush=True)
            adata = a_c
        else:
            print(f"[FILTER] WARNING: min_counts_per_cell={int(min_counts)} would drop all "
                  f"cells; skipping the counts filter.", flush=True)

    # Threshold lists start at the requested cutoff and only relax *downward*
    # (never above the user's value), so the effective per-cell gene cutoff is
    # controllable instead of a hard-coded 30.
    if min_genes is None:
        gene_thresholds = list(RELAX_GENE_THRESHOLDS)
    else:
        gene_thresholds = [v for v in sorted({int(min_genes), 20, 10, 5, 2, 1}, reverse=True)
                           if v <= int(min_genes)] or [1]
    if min_cells is None:
        cell_thresholds = list(RELAX_CELL_THRESHOLDS)
    else:
        cell_thresholds = [v for v in sorted({int(min_cells), 5, 2, 1}, reverse=True)
                           if v <= int(min_cells)] or [1]

    def genes_per_cell(a: ad.AnnData) -> np.ndarray:  # genes per cell
        return np.asarray((a.X > 0).sum(axis=1)).ravel()  # compute

    def cells_per_gene(a: ad.AnnData) -> np.ndarray:  # cells per gene
        return np.asarray((a.X > 0).sum(axis=0)).ravel()  # compute

    for mg in gene_thresholds:  # iterate gene thresholds
        gpc = genes_per_cell(adata)  # genes per cell
        keep_cells = gpc >= int(mg)  # mask
        a2 = adata[keep_cells].copy()  # subset
        if int(a2.n_obs) == 0:  # none
            continue  # next
        for mcg in cell_thresholds:  # iterate cell thresholds
            cpg = cells_per_gene(a2)  # cells per gene
            keep_genes = cpg >= int(mcg)  # mask
            a3 = a2[:, keep_genes].copy()  # subset
            if int(a3.n_obs) == 0 or int(a3.n_vars) == 0:  # invalid
                continue  # next
            print(f"[FILTER] Using MIN_GENES_PER_CELL={int(mg)}, MIN_CELLS_PER_GENE={int(mcg)} -> cells={a3.n_obs}, genes={a3.n_vars}", flush=True)  # log
            return a3  # return

    gpc0 = np.asarray((adata.X > 0).sum(axis=1)).ravel()  # original gpc
    k = int(min(int(min_cells_floor), int(adata.n_obs)))  # floor count
    top_idx = np.argsort(gpc0)[::-1][:k]  # top idx
    a_floor = adata[top_idx].copy()  # keep top
    cpg_floor = np.asarray((a_floor.X > 0).sum(axis=0)).ravel()  # gene support
    keep_genes_floor = cpg_floor > 0  # genes present
    a_floor = a_floor[:, keep_genes_floor].copy()  # subset
    print(f"[FILTER] WARNING: filtering collapsed; kept top {a_floor.n_obs} cells -> cells={a_floor.n_obs}, genes={a_floor.n_vars}", flush=True)  # warn
    return a_floor  # return


def _elbow_n_pcs(var_ratio, floor=5, cap=50):
    """Data-driven PC count: the 'knee' of the variance-ratio curve (max distance
    to the line joining first and last points). Keeps signal PCs, drops noise PCs
    that fragment the graph on targeted panels."""
    v = np.asarray(var_ratio, dtype=float)
    n = int(v.size)
    if n <= floor:
        return max(2, n)
    x = np.arange(n, dtype=float)
    x1, y1, x2, y2 = x[0], v[0], x[-1], v[-1]
    num = np.abs((y2 - y1) * x - (x2 - x1) * v + x2 * y1 - y2 * x1)
    den = np.hypot(y2 - y1, x2 - x1) + 1e-12
    k = int(np.argmax(num / den)) + 1
    return int(min(max(k, floor), cap, n))


def exact_kmeans_clustering(adata: ad.AnnData, n_clusters: int, n_pcs: int, seed: int,
                            key: str = "cluster", normalization: str = "lognorm",
                            auto_n_pcs: bool = False) -> ad.AnnData:
    X = adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)

    n_cells = int(X.shape[0])
    n_genes = int(X.shape[1])

    if n_cells == 0:
        raise ValueError("AnnData has 0 cells; cannot cluster.")

    if n_cells == 1:
        adata.obsm["X_pca_fast"] = np.zeros((1, 2), dtype=np.float32)
        adata.obs[key] = pd.Series(["0"], index=adata.obs_names, dtype=str)
        return adata

    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    Xn = (X / lib) * 1e4
    Xn = np.log1p(Xn)

    max_pcs = int(min(n_cells - 1, n_genes))
    n_pcs_eff = int(min(int(n_pcs), max_pcs))
    n_pcs_eff = int(max(2, n_pcs_eff))

    # ========================= PCA =========================
    if SCANPY_AVAILABLE:

        adata_pp = adata.copy()  # raw counts live in .X
        norm = str(normalization).lower()
        used_pearson = False

        if norm.startswith("pearson"):
            # Analytic Pearson residuals: variance-stabilizes count data, improves
            # cluster separability and marker contrast for sparse Xenium panels.
            try:
                sc.experimental.pp.normalize_pearson_residuals(adata_pp)
                used_pearson = True
                print("[PCA] Normalization: analytic Pearson residuals", flush=True)
            except Exception as e:  # older scanpy / missing experimental module
                print(f"[PCA] Pearson residuals unavailable ({e}); using lognorm.", flush=True)

        if not used_pearson:
            print("[PCA] Normalization: CP10k + log1p + scale", flush=True)
            sc.pp.normalize_total(adata_pp, target_sum=1e4)
            sc.pp.log1p(adata_pp)
            # HVG on LOG-normalized data must NOT use flavor='seurat_v3' (that flavor
            # requires raw counts and mis-selects genes here). 'seurat' is valid on
            # log data; for targeted panels this keeps ~all genes anyway.
            try:
                sc.pp.highly_variable_genes(
                    adata_pp, n_top_genes=min(2000, adata_pp.n_vars), flavor="seurat")
            except Exception:
                pass
            sc.pp.scale(adata_pp, max_value=10)

        # Number of PCs: fixed, or data-driven elbow when auto_n_pcs is on.
        cap = int(min(max(n_pcs_eff, 30) if auto_n_pcs else n_pcs_eff, max_pcs))
        cap = max(2, cap)
        sc.tl.pca(adata_pp, n_comps=cap)
        var_ratio = np.asarray(adata_pp.uns["pca"]["variance_ratio"], dtype=float)

        if auto_n_pcs:
            n_keep = _elbow_n_pcs(var_ratio, floor=5, cap=cap)
            print(f"[PCA] auto n_pcs via elbow: {n_keep} (of {cap})", flush=True)
        else:
            n_keep = min(int(n_pcs_eff), cap)

        print("[PCA] Effective PCs:", n_keep, flush=True)
        print("[PCA] Variance ratio first 10 PCs:", var_ratio[:10], flush=True)
        print("[PCA] Total variance explained (kept):",
              float(np.sum(var_ratio[:n_keep])), flush=True)

        Z = adata_pp.obsm["X_pca"][:, :n_keep]

    else:
        print("[WARNING] Scanpy not available → fallback PCA", flush=True)

        Z = PCA(n_components=n_pcs_eff, random_state=int(seed)).fit_transform(Xn)

    # ===== ALWAYS STORE PCA =====
    adata.obsm["X_pca_fast"] = Z.astype(np.float32, copy=False)

    print("[CHECK] PCA computed:", adata.obsm["X_pca_fast"].shape, flush=True)

    return adata


# Default Xenium negative-control / blank feature name patterns.
NEG_CONTROL_PATTERNS = (
    "negcontrolprobe", "negcontrolcodeword", "negprobe", "blank",
    "antisense", "unassignedcodeword", "deprecatedcodeword", "intergenic",
)


def subtract_negcontrol_background(adata: ad.AnnData,
                                   patterns=NEG_CONTROL_PATTERNS) -> ad.AnnData:
    """Estimate and subtract ambient/background signal using Xenium negative-control
    probes/codewords, then drop those control features.

    Background per cell is the mean count across control features; that expected
    background is subtracted from every real gene (floored at 0). This denoises
    profiles -> stronger, more specific markers and cleaner spatial signal. Safe
    no-op if no control features are found.
    """
    names = np.asarray([str(v).lower() for v in adata.var_names])
    is_ctrl = np.array([any(p in nm for p in patterns) for nm in names], dtype=bool)
    n_ctrl = int(is_ctrl.sum())
    if n_ctrl == 0:
        print("[BGSUB] no negative-control features found; skipping background subtraction.",
              flush=True)
        return adata

    X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=float)
    bg = X[:, is_ctrl].sum(axis=1, keepdims=True) / float(max(n_ctrl, 1))  # per-cell background
    real = ~is_ctrl
    Xr = np.clip(X[:, real] - bg, 0.0, None)

    adata2 = adata[:, real].copy()
    adata2.X = Xr.astype(np.float32)
    print(f"[BGSUB] subtracted background from {int(real.sum())} genes using "
          f"{n_ctrl} control features (mean bg/cell={float(bg.mean()):.3f}).", flush=True)
    return adata2


def remove_doublets(adata: ad.AnnData, seed: int = 0) -> ad.AnnData:
    """Flag and drop predicted doublets via Scrublet (scanpy built-in). Doublets sit
    between clusters and blunt every marker metric. Safe no-op if unavailable."""
    if not SCANPY_AVAILABLE:
        return adata
    try:
        sc.pp.scrublet(adata, random_state=int(seed))
    except Exception as e:
        print(f"[DOUBLET] Scrublet unavailable ({e}); skipping doublet removal.", flush=True)
        return adata
    if "predicted_doublet" not in adata.obs:
        return adata
    keep = ~adata.obs["predicted_doublet"].to_numpy(dtype=bool)
    n_drop = int((~keep).sum())
    if n_drop == 0 or int(keep.sum()) < 2:
        print("[DOUBLET] no doublets removed.", flush=True)
        return adata
    print(f"[DOUBLET] removed {n_drop} predicted doublets -> {int(keep.sum())} cells.",
          flush=True)
    return adata[keep].copy()
