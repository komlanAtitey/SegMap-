"""SegMap Bayesian spatial refinement and stage snapshots."""
from .common import *  # noqa: F401,F403
from .assignment import save_cell_by_gene_csv
from .clustering import silhouette_score_fast
from .markers import fast_cluster_markers_foldchange


def segmap_refinement_publication(
    adata,
    n_clusters,
    OUTPUT_DIR,
    ROI_BUDGET_FRAC=0.25,
    ROI_MAX_REGIONS=50,
    ROI_KNN=10,
    ROI_ENTROPY_PCTL=75,
    ROI_DELTA=1e-6):

        #@@@@@@@@@@ komlan import numpy as np
    import pandas as pd
    import os
    from sklearn.metrics import pairwise_distances, silhouette_samples
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix
    
    # ============================================================
    # 🔥 Initialize cluster_final (FIX)
    # ============================================================
    if "cluster_final" not in adata.obs:
        adata.obs["cluster_final"] = adata.obs["cluster"].copy()

    clusters_initial = adata.obs["cluster"].astype(int).values.copy()
    
    K = int(adata.obs["cluster"].nunique())
    Z = adata.obsm["X_pca_fast"]
    coords = adata.obs[["x_centroid_px","y_centroid_px"]].to_numpy()

    # ---------------- Degenerate case: a single cluster ----------------
    # With only one cluster there is no 2nd-nearest centroid (np.partition would
    # fail) and refinement is meaningless, so set safe defaults and return.
    if K < 2:
        print("[REFINE] WARNING: only 1 cluster found; skipping refinement. "
              "Check upstream clustering (resolution / neighbors / PCA).", flush=True)
        adata.obs["cluster_final"] = adata.obs["cluster"].astype("category")
        adata.obs["posterior_entropy"] = 0.0
        adata.obs["posterior_confidence"] = 1.0
        adata.obs["roi_selected"] = 0
        adata.obs["refinement_mode"] = "skipped_single_cluster"
        return adata

    # ---------------- Likelihood ----------------
    centroids = np.vstack([
        Z[clusters_initial==k].mean(0) if (clusters_initial==k).sum()>0 else Z.mean(0)
        for k in range(K)
    ])

    dist = pairwise_distances(Z, centroids)
    log_likelihood = -dist

    # ---------------- Spatial prior ----------------
    nn = NearestNeighbors(n_neighbors=12).fit(coords)
    neigh = nn.kneighbors(coords, return_distance=False)

    beta = min(
       1.2,
       0.3 + np.log1p(K)
    )

    # ===== REMAP CLUSTERS TO 0..K-1 =====
    adata.obs["cluster_final"] = adata.obs["cluster_final"].astype("category")
    print("[REFINE] cluster_final exists:", "cluster_final" in adata.obs, flush=True)
    
    cluster_codes = pd.Categorical(adata.obs["cluster"]).codes
    print("[REFINE] unique cluster codes:", np.unique(cluster_codes), flush=True)

    # Vectorized neighbour class-count: one-hot(K) gathered over neighbours, summed.
    onehot = np.zeros((Z.shape[0], K), dtype=np.float32)  # (n, K)
    onehot[np.arange(Z.shape[0]), cluster_codes] = 1.0  # set class membership
    log_spatial_prior = beta * onehot[neigh].sum(axis=1)  # (n, K) neighbour counts * beta

    # ---------------- Posterior ----------------
    log_posterior = log_likelihood + log_spatial_prior
    new_labels = np.argmax(log_posterior, axis=1)

    adata.obs["cluster_final"] = pd.Categorical(new_labels.astype(str))
    
    #@@@@@@@@@@ atitey
    
    # ============================================================
    # 🔥 SPATIAL SMOOTHING (CORRECT PLACE)
    # ============================================================
    from scipy.spatial import cKDTree

    coords = adata.obs[["x_centroid_px","y_centroid_px"]].to_numpy()
    tree = cKDTree(coords)

    neighbors = tree.query(coords, k=10)[1]

    cluster_final_int = adata.obs["cluster_final"].astype(int).values.copy()

    # Vectorized majority vote over k spatial neighbours (replaces per-cell loop).
    Kf = int(cluster_final_int.max()) + 1  # number of label ids present
    onehot_f = np.zeros((len(cluster_final_int), Kf), dtype=np.float32)  # (n, Kf)
    onehot_f[np.arange(len(cluster_final_int)), cluster_final_int] = 1.0  # membership
    neighbor_votes = onehot_f[neighbors].sum(axis=1)  # (n, Kf) summed neighbour votes
    cluster_final_int = neighbor_votes.argmax(axis=1).astype(int)  # majority label

    adata.obs["cluster_final"] = pd.Categorical(cluster_final_int.astype(str))
    
    #@@@@@@@@@@ atitey

    print("[REFINE] clusters updated:", adata.obs["cluster_final"].nunique(), flush=True)

    posterior = np.exp(log_posterior - log_posterior.max(1, keepdims=True))
    posterior /= posterior.sum(1, keepdims=True)

    confidence = posterior.max(1)
    entropy = -np.sum(posterior * np.log(posterior + 1e-12), 1)
    cluster_new = posterior.argmax(1)
    
    # ---------------- Common masks ----------------
    own = dist[np.arange(Z.shape[0]), clusters_initial]
    second = np.partition(dist, 1, axis=1)[:, 1]
    margin = second - own
    boundary_mask = (margin < np.percentile(margin, 35))

    # NOTE: a per-cell silhouette pass used to be computed here but its result
    # (sil_mask) was never consumed; removed to save an O(n^2) computation.

    improvement = (
        posterior[np.arange(len(cluster_new)), cluster_new] -
        posterior[np.arange(len(cluster_new)), clusters_initial]
    )

    # ===== RELAXED BASE MASK (FIXED) =====
    base_mask = (
        (confidence > 0.6) &  #@@@@@@@@@ komlan 0.5) &
        (improvement > 0.0) & #@@@@@@@@@ komlan 0.02) &
        (cluster_new != clusters_initial)
    )
    
    # ============================================================
    # 🔹 OPTION 1: GLOBAL REFINEMENT (original behavior)
    # ============================================================
    cluster_global = clusters_initial.copy()
    mask_global = base_mask.copy()
    max_updates = int(0.5 * len(cluster_global)) #@@@@@@@@ 0.6
    idx = np.where(mask_global)[0]

    if len(idx) > max_updates:
        keep = idx[np.argsort(confidence[idx])[-max_updates:]]
        mask_global[:] = False
        mask_global[keep] = True

    cluster_global[mask_global] = cluster_new[mask_global]

    # ============================================================
    # 🔹 OPTION 2: ROI REFINEMENT
    # ============================================================

    # ROI detection
    entropy_thresh = np.percentile(entropy, ROI_ENTROPY_PCTL)
    high_uncertainty = entropy >= entropy_thresh

    nn_roi = NearestNeighbors(n_neighbors=ROI_KNN).fit(coords)
    neigh_roi = nn_roi.kneighbors(coords, return_distance=False)

    rows, cols = [], []
    for i in range(len(coords)):
        if not high_uncertainty[i]:
            continue
        for j in neigh_roi[i]:
            if high_uncertainty[j]:
                rows.append(i)
                cols.append(j)

    graph = csr_matrix((np.ones(len(rows)), (rows, cols)),
                       shape=(len(coords), len(coords)))

    n_comp, labels = connected_components(graph, directed=False)

    roi_list = []
    for cid in range(n_comp):
        idx = np.where(labels == cid)[0]
        if len(idx) == 0:
            continue
        if not np.any(high_uncertainty[idx]):
            continue

        utility = entropy[idx].mean()
        cost = len(idx)
        rho = utility / max(cost, ROI_DELTA)

        roi_list.append((rho, idx))

    roi_list.sort(key=lambda t: t[0], reverse=True)

    budget = int(ROI_BUDGET_FRAC * len(coords))
    selected_cells = np.zeros(len(coords), dtype=bool)

    used = 0
    for rho, idx in roi_list:
        if used + len(idx) > budget:
            continue
        selected_cells[idx] = True
        used += len(idx)
        if used >= budget:
            break

    # ROI mask
    mask_roi = base_mask & selected_cells
    
    # ===== DEBUG (ADD HERE) =====
    print("Base mask count:", base_mask.sum(), flush=True)
    print("ROI mask count:", mask_roi.sum(), flush=True)

    cluster_roi = clusters_initial.copy()
    cluster_roi[mask_roi] = cluster_new[mask_roi]

    # ============================================================
    # AUTO-SELECTION
    # ============================================================
    def score(clusters):
        sil_mean = silhouette_score_fast(Z, clusters)  # subsampled, fast
        sil_mean = 0.0 if not np.isfinite(sil_mean) else sil_mean
        return (
            -entropy.mean() +      # lower entropy is better
            sil_mean +             # better structure
            confidence.mean()      # higher certainty
        )
    # ------------------------------------------------------------
    # Compare ROI vs GLOBAL
    # ------------------------------------------------------------
    score_global = score(cluster_global)
    score_roi = score(cluster_roi)

    print(
        f"[REFINEMENT] GLOBAL score = {score_global:.4f}",
        flush=True
    )
    print(
        f"[REFINEMENT] ROI score = {score_roi:.4f}",
        flush=True
    )
    if score_roi > score_global:
        cluster_final = cluster_roi
        selected_mode = "ROI"
        final_mask = mask_roi
    else:
        cluster_final = cluster_global
        selected_mode = "GLOBAL"
        final_mask = mask_global


    # ============================================================
    # DIAGNOSTIC: Compare clustering before vs after refinement
    # ============================================================
    try:
        score_initial = silhouette_score_fast(Z, clusters_initial)
        score_final = silhouette_score_fast(Z, cluster_final)
        print(
            f"[REFINEMENT] Initial silhouette = {score_initial:.4f}",
            flush=True
        )

        print(
            f"[REFINEMENT] Final silhouette = {score_final:.4f}",
            flush=True
        )
    except Exception as e:
        print(
            f"[REFINEMENT] Silhouette comparison failed: {e}",
            flush=True
        )


    # ============================================================
    # SAVE RESULTS
    # ============================================================
    adata.obs["cluster_final"] = pd.Categorical(
        cluster_final.astype(str)
    )
    adata.obs["posterior_entropy"] = entropy
    adata.obs["posterior_confidence"] = confidence
    adata.obs["roi_selected"] = selected_cells.astype(int)
    adata.obs["refinement_mode"] = selected_mode

    # ===== FIX 5: DEBUG HOW MANY CELLS CHANGED =====
    #@@@@@@@@
    # ============================================================
    # Refinement diagnostics
    # ============================================================
    changed = int((clusters_initial != cluster_final).sum())

    n_cells = int(len(cluster_final))

    pct_changed = 100.0 * changed / max(1, n_cells)

    pct_retained = 100.0 - pct_changed

    print(
        f"[REFINEMENT] Cells changed: "
        f"{changed}/{n_cells} "
        f"({pct_changed:.2f}%)",
        flush=True
    )

    print(
        f"[REFINEMENT] Cells retained: "
        f"{n_cells - changed}/{n_cells} "
        f"({pct_retained:.2f}%)",
        flush=True
    )

    # ============================================================
    # Safety check
    # ============================================================
    MAX_CHANGED_FRAC = 0.30   # allow at most 30% reassignment
    frac_changed = changed / max(1, n_cells)
    if frac_changed > MAX_CHANGED_FRAC:
        print(
            "[REFINEMENT] WARNING: excessive reassignment "
            f"detected ({pct_changed:.2f}%). "
            "Reverting to initial clustering.",
            flush=True
        )
        cluster_final = clusters_initial.copy()
        changed = 0
        print(
            "[REFINEMENT] Reverted to initial clustering.",
            flush=True
        )
    else:
        print(
            "[REFINEMENT] Refinement accepted.",
            flush=True
        )
    #@@@@@@@@
    print(f"[REFINEMENT] Cells changed: {changed}", flush=True)

    pd.DataFrame({
        "cell_id": adata.obs_names,
        "entropy": entropy,
        "confidence": confidence,
        "in_roi": selected_cells.astype(int),
        "updated": final_mask.astype(int)
    }).to_csv(
        os.path.join(OUTPUT_DIR,"posterior_entropy_final.csv"),
        index=False
    )

    return adata


def save_segmap_stage(adata,FIG_DIR,OUTPUT_DIR,suffix):

    import os
    import pandas as pd

    markers=fast_cluster_markers_foldchange(
        adata,"cluster",10
    )

    markers.to_csv(
        os.path.join(
            FIG_DIR,
            f"cluster_markers_all_fast_{suffix}.csv"
        ),
        index=False
    )

    spatial_df=adata.obs[[
        "cluster",
        "x_centroid_px",
        "y_centroid_px"
    ]].copy()

    spatial_df.insert(
        0,"cell_id",adata.obs_names
    )

    spatial_df.to_csv(
        os.path.join(
            FIG_DIR,
            f"spatial_clusters_{suffix}.csv"
        ),
        index=False
    )

    if "X_umap_fast" in adata.obsm:

        Y=adata.obsm["X_umap_fast"]

        pd.DataFrame({
            "cell_id":adata.obs_names,
            "umap1":Y[:,0],
            "umap2":Y[:,1],
            "cluster":adata.obs["cluster"]
        }).to_csv(
            os.path.join(
                FIG_DIR,
                f"umap_clusters_{suffix}.csv"
            ),
            index=False
        )

    save_cell_by_gene_csv(
        adata,
        os.path.join(
            OUTPUT_DIR,
            f"cell_by_gene_for_umap_{suffix}.csv"
        )
    )
