"""Silhouette helpers and Leiden clustering with resolution tuning."""
from .common import *  # noqa: F401,F403


def silhouette_score_fast(X, labels, max_cells: int = 10000, seed: int = 0) -> float:
    """Silhouette score on a random subsample.

    silhouette_score is O(n^2) in both time and memory; on large Xenium runs it
    was being recomputed on every cell several times. Subsampling to `max_cells`
    keeps the estimate stable while turning minutes of compute into seconds.
    """
    from sklearn.metrics import silhouette_score as _sil
    X = np.asarray(X)
    labels = np.asarray(labels)
    n = X.shape[0]
    if n > int(max_cells):
        rng = np.random.default_rng(int(seed))
        sel = rng.choice(n, size=int(max_cells), replace=False)
        X, labels = X[sel], labels[sel]
    # silhouette needs >=2 labels and each label present; guard against degenerate subsamples
    uniq = np.unique(labels)
    if uniq.size < 2:
        return float("nan")
    return float(_sil(X, labels))


def silhouette_samples_fast(X, labels, max_cells: int = 10000, seed: int = 0):
    """Per-cell silhouette, but computed on a subsample for large datasets.

    Returns a full-length float array; cells outside the subsample are NaN
    (use np.nanmean / np.nanmedian on the result). For n <= max_cells this is
    the exact per-cell silhouette.
    """
    from sklearn.metrics import silhouette_samples as _sils
    X = np.asarray(X)
    labels = np.asarray(labels)
    n = X.shape[0]
    out = np.full(n, np.nan, dtype=np.float32)
    if np.unique(labels).size < 2:
        return out
    if n > int(max_cells):
        rng = np.random.default_rng(int(seed))
        sel = rng.choice(n, size=int(max_cells), replace=False)
        sub = silhouette_samples_fast(X[sel], labels[sel], max_cells=max_cells, seed=seed)
        out[sel] = sub
        return out
    out[:] = _sils(X, labels).astype(np.float32)
    return out


def run_leiden_best(adata, base_resolution, tune=True, key="cluster_leiden",
                    seed=0, rep="X_pca_fast", min_clusters=2):
    """Run Leiden with the fast igraph backend and (optionally) tune resolution.

    Cluster-separation upgrades:
      * Uses flavor="igraph", n_iterations=2, directed=False (the modern Scanpy
        recommendation): faster and yields better-modularity partitions than the
        legacy leidenalg default.
      * When tune=True, sweeps resolutions around the requested value and keeps the
        partition with the highest (subsampled) silhouette, so clusters are chosen
        for separation rather than an arbitrary fixed resolution.
    Requires an existing neighbours graph (sc.pp.neighbors) and adata.obsm[rep].
    """
    def _leiden(res, key_added, force_legacy=False):
        if not force_legacy:
            try:
                sc.tl.leiden(adata, resolution=res, key_added=key_added,
                             flavor="igraph", n_iterations=2, directed=False,
                             random_state=int(seed))
                return
            except TypeError:  # older scanpy without igraph flavor kwargs
                pass
        # legacy leidenalg backend (used as a fallback when igraph collapses)
        sc.tl.leiden(adata, resolution=res, key_added=key_added,
                     random_state=int(seed))

    def _ncl(keyname):
        return int(adata.obs[keyname].astype(str).nunique())

    def _escalate_to_min(res, key_added, min_clusters):
        """Run Leiden; if it yields fewer than `min_clusters`, raise the resolution
        (and finally try the legacy backend) until it does, tracking the best result.
        A too-low resolution on a well-connected graph collapses all structure into
        1-2 clusters even when the UMAP clearly shows many groups; this prevents that."""
        _leiden(res, key_added)
        n = _ncl(key_added)
        if n >= int(min_clusters):
            return res
        best_n, best_res = n, res
        for mult in (1.5, 2.0, 3.0, 4.0, 6.0, 8.0):  # climb resolution
            r = res * mult
            _leiden(r, key_added)
            n = _ncl(key_added)
            if n > best_n:
                best_n, best_res = n, r
            if n >= int(min_clusters):
                print(f"[CLUSTER] resolution escalated {res:g} -> {r:g} "
                      f"({n} clusters; target >= {int(min_clusters)}).", flush=True)
                return r
        _leiden(best_res, key_added, force_legacy=True)  # try legacy backend
        if _ncl(key_added) > best_n:
            print(f"[CLUSTER] used leidenalg at resolution={best_res:g} -> "
                  f"{_ncl(key_added)} clusters.", flush=True)
            return best_res
        _leiden(best_res, key_added)  # restore best igraph partition
        print(f"[CLUSTER] WARNING: reached only {best_n} clusters (target >= "
              f"{int(min_clusters)}); raise --resolution for finer clustering.", flush=True)
        return best_res

    if not tune:
        _escalate_to_min(base_resolution, key, min_clusters)
        return float(base_resolution)

    Z = adata.obsm.get(rep)
    # include upward multipliers so a too-low base can still reach real structure
    candidates = sorted(set(round(float(base_resolution) * m, 3)
                            for m in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)))
    results = []  # (res, k, silhouette, labels)
    for res in candidates:
        _leiden(res, "_leiden_tmp")
        labels = adata.obs["_leiden_tmp"].astype(int).to_numpy()
        k = int(np.unique(labels).size)
        if k < 2:
            print(f"[CLUSTER] resolution={res} -> {k} cluster (skipped)", flush=True)
            continue
        s = silhouette_score_fast(Z, labels, seed=int(seed)) if Z is not None else np.nan
        print(f"[CLUSTER] resolution={res} -> {k} clusters, silhouette={s:.4f}", flush=True)
        results.append((res, k, s, labels.copy()))

    if "_leiden_tmp" in adata.obs:
        del adata.obs["_leiden_tmp"]
    if not results:  # everything collapsed; escalate + fall back
        _escalate_to_min(base_resolution, key, min_clusters)
        return float(base_resolution)

    # Silhouette alone is biased toward very few clusters, so gate on count first:
    # prefer partitions with >= min_clusters; among those pick the best silhouette.
    viable = [r for r in results if r[1] >= int(min_clusters)]
    if viable:
        best_res, best_k, best_s, best_labels = max(
            viable, key=lambda r: (r[2] if np.isfinite(r[2]) else -np.inf))
    else:  # none reached the floor -> take the most-clustered partition
        best_res, best_k, best_s, best_labels = max(results, key=lambda r: r[1])
    adata.obs[key] = pd.Categorical(best_labels.astype(str))
    print(f"[CLUSTER] Selected resolution={best_res} "
          f"(silhouette={best_s:.4f}, k={best_k})", flush=True)
    return float(best_res)


def spatial_denoise_labels(labels, coords, confidence=None, k=15,
                           agree_frac=0.65, conf_pctl=50, n_iter=2, seed=0):
    """Confidence-gated spatial denoising of categorical labels.

    Keeps transcriptomic cell identity intact and only cleans up isolated
    misassignments: a cell is relabeled to its spatial-neighborhood majority
    *only* when (a) the cell is low-confidence (posterior confidence at/below the
    `conf_pctl` percentile; if `confidence` is None, all cells are eligible) AND
    (b) at least `agree_frac` of its `k` nearest spatial neighbors share a single
    label that differs from the cell's own. Confident / interior cells are never
    touched, so genuine boundaries and identity are preserved while salt-and-
    pepper noise is removed. Returns an integer label array.
    """
    labels = np.asarray(labels).astype(int).copy()
    coords = np.asarray(coords, dtype=float)
    n = labels.shape[0]
    if int(n_iter) <= 0 or n < 3:
        return labels

    from scipy.spatial import cKDTree
    kq = int(min(k + 1, n))
    nbr = cKDTree(coords).query(coords, k=kq)[1]
    if nbr.ndim == 1:
        nbr = nbr[:, None]
    if nbr.shape[1] >= 2:
        nbr = nbr[:, 1:]  # drop self (nearest neighbor is the cell itself)

    if confidence is not None:
        confidence = np.asarray(confidence, dtype=float)
        thr = np.percentile(confidence, conf_pctl)
        eligible = confidence <= thr
    else:
        eligible = np.ones(n, dtype=bool)

    kk = float(nbr.shape[1])
    for _ in range(int(n_iter)):
        K = int(labels.max()) + 1
        onehot = np.zeros((n, K), dtype=np.float32)
        onehot[np.arange(n), labels] = 1.0
        votes = onehot[nbr].sum(axis=1)            # (n, K) neighbor label counts
        maj = votes.argmax(axis=1)                 # majority neighbor label
        maj_frac = votes.max(axis=1) / kk          # its agreement fraction
        change = eligible & (maj != labels) & (maj_frac >= float(agree_frac))
        if not np.any(change):
            break
        labels[change] = maj[change]
    return labels


def cluster_separability_report(adata, cluster_key="cluster_final",
                                rep="X_pca_fast", out_csv=None, seed=0):
    """Quantify how separated each cluster is, so 'overlap' can be interpreted.

    For each cluster reports: mean silhouette in `rep` space (low/negative =
    mixed with others), its nearest cluster by centroid distance, and that
    distance. Clusters with low silhouette whose nearest neighbour shares top
    markers are candidates to MERGE (over-split); a smooth gradient of moderate
    silhouettes is a real continuum; an isolated low-silhouette cluster is likely
    low-quality/doublets. Writes a CSV and prints a summary.
    """
    if rep not in adata.obsm:
        return None
    Z = np.asarray(adata.obsm[rep])
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    codes = pd.Categorical(labels).codes
    sil = silhouette_samples_fast(Z, codes, seed=int(seed))

    names = sorted(set(labels), key=lambda s: (len(s), s))
    cents = np.vstack([Z[labels == c].mean(axis=0) for c in names])
    D = np.sqrt(((cents[:, None, :] - cents[None, :, :]) ** 2).sum(-1))

    rows = []
    for i, c in enumerate(names):
        m = labels == c
        s = float(np.nanmean(sil[m])) if m.any() else float("nan")
        d = D[i].copy(); d[i] = np.inf
        j = int(np.argmin(d))
        rows.append({"cluster": c, "n": int(m.sum()), "mean_silhouette": s,
                     "nearest_cluster": names[j], "centroid_dist": float(D[i, j])})
    df = pd.DataFrame(rows).sort_values("mean_silhouette").reset_index(drop=True)
    if out_csv:
        df.to_csv(out_csv, index=False)
    print("[SEPARABILITY] per-cluster (low/negative silhouette = overlapping):", flush=True)
    for _, r in df.iterrows():
        print(f"   cluster {str(r['cluster']):<4} n={int(r['n']):<6} "
              f"sil={r['mean_silhouette']:+.3f}  nearest={r['nearest_cluster']} "
              f"(d={r['centroid_dist']:.2f})", flush=True)
    return df


def merge_indistinct_clusters(labels, Z, threshold=1.0, min_clusters=2, max_iter=200):
    """Consolidate over-split clusters that are not separable in `Z` (PCA space).

    Separation of a pair is the ratio of between-centroid distance to the sum of
    within-cluster radii; pairs below `threshold` overlap and are merged (closest
    first), iterating until every remaining pair is separated or `min_clusters`
    is reached. Fewer, larger, purer clusters -> more cells per cluster -> stronger
    DE significance, while improving interpretability. Cheap (centroids + radii).
    Returns integer labels relabelled 0..K-1.
    """
    labels = np.asarray(labels).astype(int).copy()
    Z = np.asarray(Z, dtype=float)
    for _ in range(int(max_iter)):
        uniq = np.unique(labels)
        if len(uniq) <= int(min_clusters):
            break
        cents, radii = {}, {}
        for c in uniq:
            pts = Z[labels == c]
            cents[c] = pts.mean(axis=0)
            radii[c] = float(np.linalg.norm(pts - cents[c], axis=1).mean()) if len(pts) > 1 else 0.0
        best = None  # (sep, a, b)
        for i in range(len(uniq)):
            for j in range(i + 1, len(uniq)):
                a, b = uniq[i], uniq[j]
                d = float(np.linalg.norm(cents[a] - cents[b]))
                sep = d / (radii[a] + radii[b] + 1e-9)
                if best is None or sep < best[0]:
                    best = (sep, a, b)
        if best is None or best[0] >= float(threshold):
            break
        _, a, b = best
        labels[labels == b] = a  # merge smaller-id target
    uniq = np.unique(labels)
    remap = {c: i for i, c in enumerate(uniq)}
    return np.array([remap[v] for v in labels], dtype=int)
