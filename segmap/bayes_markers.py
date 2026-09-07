"""
bayes_markers.py — Bayesian marker model for SegMap (Fig. 1D).

Implements the per-gene spatially-variable-gene model by Bayesian model selection
between a spatial (Gaussian-process / graph-spectral) model and a non-spatial null.
This is the *implemented* counterpart of the specification: it produces, per gene,
a posterior marker probability p_j and an effect size FSV_j, plus cluster-specific
markers (spatial posterior gated by one-vs-rest enrichment).

Scalable graph-spectral variant: the smooth low-frequency eigenvectors of the kNN
graph Laplacian act as the spatial basis, so each gene is fit in O(N*m) after a
one-time sparse partial eigendecomposition (m << N). Reuses SegMap's spatial graph.

Dependencies: numpy, scipy, scikit-learn (already required by the pipeline).
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, diags, identity
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import NearestNeighbors

__all__ = [
    "build_spatial_basis", "fit_all_genes", "select_markers_bfdr",
    "cluster_bayes_markers", "run_bayes_markers",
]


# --------------------------------------------------------------------------- #
# 1. Spatial basis: low-frequency Laplacian eigenmaps of the kNN graph
# --------------------------------------------------------------------------- #
def build_spatial_basis(coords: np.ndarray, k: int = 15, m: int = 100,
                        normalized: bool = True):
    """Return the m smoothest non-trivial Laplacian eigenvectors U (N x m) and
    eigenvalues gamma (m,) of a symmetric kNN graph on `coords`."""
    coords = np.asarray(coords, dtype=float)
    N = coords.shape[0]
    m = int(min(m, max(1, N - 2)))
    k = int(min(k, max(1, N - 1)))

    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords)
    idx = idx[:, 1:]  # drop self
    rows = np.repeat(np.arange(N), k)
    cols = idx.reshape(-1)
    data = np.ones(rows.size)
    W = csr_matrix((data, (rows, cols)), shape=(N, N))
    W = ((W + W.T) > 0).astype(float)          # symmetrize -> 0/1 adjacency
    deg = np.asarray(W.sum(1)).ravel()

    if normalized:
        dinv = 1.0 / np.sqrt(np.maximum(deg, 1e-8))
        Dinv = diags(dinv)
        L = identity(N) - Dinv @ W @ Dinv       # normalized Laplacian (PSD, spectrum [0,2])
        sigma = -1e-6
    else:
        L = diags(deg) - W                       # combinatorial Laplacian
        sigma = -1e-6

    # smallest eigenpairs; shift-invert around ~0 for the low-frequency modes
    vals, vecs = eigsh(L.tocsc(), k=min(m + 1, N - 1), sigma=sigma, which="LM")
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    keep = np.where(vals > 1e-8)[0][:m]          # drop trivial ~0 mode(s)
    return {"U": vecs[:, keep], "gamma": vals[keep], "N": N, "m": len(keep)}


# --------------------------------------------------------------------------- #
# 2. Per-gene evidence under spatial vs null model
# --------------------------------------------------------------------------- #
def _nll_lowrank(par, yt2, gamma, ss_perp, N, m):
    s2, e2 = np.exp(par)                          # positive via log-parameterization
    d = np.maximum(s2 / gamma + e2, 1e-10)        # variance of the m projections
    return 0.5 * (np.sum(yt2 / d + np.log(d)) + ss_perp / e2
                  + (N - m) * np.log(e2) + N * np.log(2 * np.pi))


def _fit_gene(y, basis):
    U, gamma, N, m = basis["U"], basis["gamma"], basis["N"], basis["m"]
    yc = y - y.mean()
    yt = U.T @ yc
    yt2 = yt * yt
    ss_perp = max(float(yc @ yc) - float(yt2.sum()), 0.0)
    v = max(float(yc @ yc) / N, 1e-8)
    opt = minimize(_nll_lowrank, x0=np.log([v / 2, v / 2]),
                   args=(yt2, gamma, ss_perp, N, m), method="L-BFGS-B")
    s2, e2 = np.exp(opt.x)
    logL1 = -opt.fun
    logL0 = -0.5 * N * (np.log(2 * np.pi * v) + 1.0)   # null: s2 = 0
    return logL1, logL0, s2, e2, s2 / (s2 + e2)


def fit_all_genes(expr: np.ndarray, basis, dpar: int = 1, pi0: float = 0.5,
                  eb_iter: int = 100, eb_tol: float = 1e-4):
    """Fit every gene (expr: cells x genes). Returns dict of arrays:
    logBF, posterior (p_j), fsv, sigma_s2, sigma_e2, pi_hat."""
    expr = np.asarray(expr, dtype=float)
    G = expr.shape[1]
    N = basis["N"]
    logL1 = np.empty(G); logL0 = np.empty(G)
    fsv = np.empty(G); s2 = np.empty(G); e2 = np.empty(G)
    for j in range(G):
        logL1[j], logL0[j], s2[j], e2[j], fsv[j] = _fit_gene(expr[:, j], basis)
    logBF = logL1 - logL0 - 0.5 * dpar * np.log(N)     # BIC-penalized log Bayes factor

    def sigmoid(x):
        return 1.0 / (1.0 + np.exp(-x))
    logit = lambda p: np.log(p / (1 - p))

    pi_hat = pi0                                        # empirical-Bayes prior
    for _ in range(eb_iter):
        p = sigmoid(logBF + logit(pi_hat))
        new = min(max(float(p.mean()), 1e-4), 1 - 1e-4)
        if abs(new - pi_hat) < eb_tol:
            pi_hat = new; break
        pi_hat = new
    p = sigmoid(logBF + logit(pi_hat))
    return {"logBF": logBF, "posterior": p, "fsv": fsv,
            "sigma_s2": s2, "sigma_e2": e2, "pi_hat": pi_hat}


def select_markers_bfdr(posterior: np.ndarray, alpha: float = 0.05):
    """Largest set with mean(1 - p) <= alpha (Bayesian FDR). Returns index array."""
    posterior = np.asarray(posterior)
    order = np.argsort(-posterior)
    fdr = np.cumsum(1.0 - posterior[order]) / (np.arange(order.size) + 1)
    keep = np.where(fdr <= alpha)[0]
    return order[: keep.max() + 1] if keep.size else np.array([], dtype=int)


# --------------------------------------------------------------------------- #
# 3. Cluster-specific Bayesian markers (spatial posterior gated by enrichment)
# --------------------------------------------------------------------------- #
def cluster_bayes_markers(expr, fit, clusters, genes, alpha=0.05, min_contrast=0.0):
    """Return {cluster: [marker genes]} — genes that are spatial markers (Bayesian
    FDR) AND positively enriched one-vs-rest in the cluster."""
    expr = np.asarray(expr, dtype=float)
    clusters = np.asarray(clusters)
    sel = select_markers_bfdr(fit["posterior"], alpha)
    sel_mask = np.zeros(expr.shape[1], dtype=bool); sel_mask[sel] = True
    out = {}
    for c in np.unique(clusters):
        inC = clusters == c
        if inC.sum() == 0 or (~inC).sum() == 0:
            out[str(c)] = []; continue
        contrast = expr[inC].mean(0) - expr[~inC].mean(0)
        keep = np.where(sel_mask & (contrast > min_contrast))[0]
        keep = keep[np.argsort(-fit["posterior"][keep])]
        out[str(c)] = [genes[i] for i in keep]
    return out


# --------------------------------------------------------------------------- #
# 4. Top-level entry: run on an AnnData and write outputs
# --------------------------------------------------------------------------- #
def run_bayes_markers(adata, cluster_key="cluster_final", output_dir=".",
                      coord_keys=("x_centroid_px", "y_centroid_px"),
                      k=15, m=100, alpha=0.05, max_cells=None, seed=0):
    """Fit the Bayesian marker model on `adata` and write:
        Benchmark_bayes_markers.csv        (per-gene p_j, FSV, logBF)
        Bayes_cluster_markers.csv          (per-cluster marker genes)
    Returns the per-gene results dict. Safe to call inside the pipeline; on any
    failure it logs and returns None rather than aborting the run.
    """
    import os
    try:
        import pandas as pd
        rng = np.random.default_rng(seed)

        X = adata.X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        X = X.astype(float)
        coords = adata.obs[list(coord_keys)].to_numpy(dtype=float)
        genes = list(adata.var_names)
        clusters = (adata.obs[cluster_key].to_numpy()
                    if cluster_key in adata.obs else None)

        # optional subsample for the (one-time) eigendecomposition on huge sections
        if max_cells is not None and X.shape[0] > max_cells:
            idx = rng.choice(X.shape[0], size=int(max_cells), replace=False)
            X, coords = X[idx], coords[idx]
            clusters = clusters[idx] if clusters is not None else None

        basis = build_spatial_basis(coords, k=k, m=m)
        fit = fit_all_genes(X, basis)

        df = pd.DataFrame({
            "gene": genes,
            "posterior_marker_prob": fit["posterior"],
            "fsv": fit["fsv"],
            "logBF": fit["logBF"],
            "sigma_s2": fit["sigma_s2"],
            "sigma_e2": fit["sigma_e2"],
        }).sort_values("posterior_marker_prob", ascending=False)
        df.attrs["pi_hat"] = fit["pi_hat"]
        path1 = os.path.join(output_dir, "Benchmark_bayes_markers.csv")
        df.to_csv(path1, index=False)
        print(f"[BAYES] wrote {path1}  (pi_hat={fit['pi_hat']:.3f}, "
              f"n_markers@FDR{alpha}={select_markers_bfdr(fit['posterior'], alpha).size})",
              flush=True)

        if clusters is not None:
            cm = cluster_bayes_markers(X, fit, clusters, genes, alpha=alpha)
            rows = [{"cluster": c, "rank": r + 1, "gene": g}
                    for c, gs in cm.items() for r, g in enumerate(gs)]
            path2 = os.path.join(output_dir, "Bayes_cluster_markers.csv")
            pd.DataFrame(rows).to_csv(path2, index=False)
            print(f"[BAYES] wrote {path2}", flush=True)

        return fit
    except Exception as e:  # never abort the main pipeline on the marker add-on
        print(f"[BAYES] skipped (non-fatal): {type(e).__name__}: {e}", flush=True)
        return None


# --------------------------------------------------------------------------- #
# 5. Close the D->E loop: couple gene-relevance p_j into the cell-state
#    posterior q_i(phi). Genes enter the label feature space weighted by p_j,
#    so q_i(phi)=p(z_i=phi|X,s) and its entropy H_i inherit marker uncertainty.
# --------------------------------------------------------------------------- #
def marker_weighted_pca(X, posterior, n_pcs=50):
    """PCA of expression with each gene scaled by sqrt(p_j) (soft marker gating)."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    X = np.asarray(X, dtype=float)
    w = np.sqrt(np.clip(np.asarray(posterior, float), 0.0, 1.0))
    Xw = StandardScaler(with_std=True).fit_transform(X) * w[None, :]
    n_pcs = int(min(n_pcs, min(Xw.shape) - 1))
    return PCA(n_components=n_pcs, random_state=0).fit_transform(Xw)


def _potts_log_prior(labels, coords, K, beta, k=12):
    """Potts prior contribution beta * (#neighbours in class k), over a kNN graph."""
    labels = np.asarray(labels)
    N = labels.shape[0]
    k = int(min(k, max(1, N - 1)))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    _, idx = nn.kneighbors(coords); idx = idx[:, 1:]
    lp = np.zeros((N, K))
    for c in range(K):
        lp[:, c] = beta * (labels[idx] == c).sum(1)
    return lp


def couple_pj_to_qi(adata, fit, cluster_key="cluster_final", n_pcs=50, k=12,
                    coord_keys=("x_centroid_px", "y_centroid_px")):
    """Recompute the cell-state posterior q_i(phi) and entropy H_i on a
    marker-weighted feature space (conditioned on X via p_j), closing the
    D->E loop. Writes obs['posterior_entropy_coupled'] and
    obs['posterior_confidence_coupled']; returns (q, H, confidence)."""
    import pandas as pd
    X = adata.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    coords = adata.obs[list(coord_keys)].to_numpy(dtype=float)
    labels = pd.Categorical(adata.obs[cluster_key]).codes
    K = int(labels.max() + 1)
    q, H, conf = _coupled_posterior_arrays(X, fit["posterior"], labels, coords,
                                           n_pcs=n_pcs, k=k)
    adata.obs["posterior_entropy_coupled"] = H
    adata.obs["posterior_confidence_coupled"] = conf
    return q, H, conf


def _coupled_posterior_arrays(X, posterior, labels, coords, n_pcs=50, k=12):
    """Core array implementation (also used for standalone validation)."""
    from sklearn.metrics import pairwise_distances
    Zw = marker_weighted_pca(X, posterior, n_pcs=n_pcs)
    K = int(labels.max() + 1)
    centroids = np.vstack([Zw[labels == c].mean(0) if (labels == c).any()
                           else Zw.mean(0) for c in range(K)])
    log_lik = -pairwise_distances(Zw, centroids)              # Gaussian-like likelihood
    beta = min(1.2, 0.3 + np.log(1 + K))                     # same schedule as refinement
    log_prior = _potts_log_prior(labels, coords, K, beta, k=k)
    log_post = log_lik + log_prior
    q = np.exp(log_post - log_post.max(1, keepdims=True))
    q /= q.sum(1, keepdims=True)
    conf = q.max(1)
    H = -np.sum(q * np.log(q + 1e-12), 1)                    # cell-state entropy
    return q, H, conf
