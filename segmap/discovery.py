"""
discovery.py — SegMap uncertainty-guided DIAGNOSTIC layer (Fig. 1d).

Read-only diagnosis and prioritization. This module propagates marker-inference
uncertainty—via the implemented forward coupling p_j -> q_i
(bayes_markers.couple_pj_to_qi)—into per-cell annotation confidence, and surfaces
biologically informative regions: ambiguous cell states, transitional phenotypes,
tissue boundaries and candidate cell-cell interaction niches. It is DIAGNOSTIC,
not corrective: it does NOT re-estimate markers, re-assign cells, alter cell
definitions, or feed back into reconstruction. The iterative
p_j -> q_i -> p_j loop (re-extracting transcripts within high-uncertainty ROIs and
re-estimating cell identities) is a planned extension and is not implemented here.
Outputs are written to a separate folder and never overwrite the reconstruction or
marker outputs the reported results depend on; the module is off by default.

Metrics follow Fig. 1d exactly:
  entropy         H_i = -sum_phi q_i(phi) log q_i(phi)
  boundary score  M_i = 1 - max_phi q_i(phi)
  interaction     S_ij = (w_ij / d_ij^alpha) * s_ij   (s_ij = cosine transcriptomic similarity)
  ROI utility     rho_i = H_i / max(C_i, delta)       (C_i = transcripts in ROI Q_i)
Overlapping ROIs are RETAINED (not merged) and each is characterized in place.
"""
from __future__ import annotations
import os, json
import numpy as np
from dataclasses import dataclass, field
from scipy.sparse import coo_matrix, csr_matrix
from sklearn.neighbors import NearestNeighbors

__all__ = ["ROI", "cell_state_entropy", "boundary_score", "build_spatial_weights",
           "interaction_score", "build_rois", "characterize_roi", "run_discovery"]


@dataclass
class ROI:
    center_cell: int
    bounds: tuple            # (x0, x1, y0, y1) in the same units as centroids
    entropy: float
    n_transcripts: int
    utility_cost: float
    member_cells: np.ndarray = field(default_factory=lambda: np.array([], int))


def cell_state_entropy(q: np.ndarray) -> np.ndarray:
    """Per-cell Shannon entropy H_i of the categorical cell-state posterior q
    (cells x states)."""
    q = np.clip(np.asarray(q, float), 1e-12, 1.0)
    q = q / q.sum(1, keepdims=True)
    return -np.sum(q * np.log(q), axis=1)


def boundary_score(q: np.ndarray) -> np.ndarray:
    """Boundary/transition score M_i = 1 - max_phi q_i(phi)."""
    q = np.clip(np.asarray(q, float), 0.0, 1.0)
    return 1.0 - q.max(axis=1)


def build_spatial_weights(centroids: np.ndarray, k: int = 15) -> csr_matrix:
    """Symmetric kNN adjacency over centroids (used only if a graph is not
    already available)."""
    N = centroids.shape[0]
    k = int(min(k, max(1, N - 1)))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(centroids)
    _, idx = nn.kneighbors(centroids)
    idx = idx[:, 1:]
    rows = np.repeat(np.arange(N), k)
    cols = idx.reshape(-1)
    A = coo_matrix((np.ones(rows.size), (rows, cols)), shape=(N, N)).tocsr()
    return ((A + A.T) > 0).astype(float)


def interaction_score(W: csr_matrix, centroids: np.ndarray, expr: np.ndarray,
                      alpha: float = 1.0, eps: float = 1e-8) -> csr_matrix:
    """Candidate interaction-niche score over the edges of W (Fig. 1d):
    S_ij = (w_ij / d_ij^alpha) * s_ij, with s_ij the non-negative cosine
    transcriptomic similarity."""
    Wc = W.tocoo()
    i, j = Wc.row, Wc.col
    d = np.linalg.norm(centroids[i] - centroids[j], axis=1) + eps
    X = np.asarray(expr, float)
    Xn = X / (np.linalg.norm(X, axis=1, keepdims=True) + eps)
    s = np.einsum("ij,ij->i", Xn[i], Xn[j])
    S = (Wc.data / np.power(d, alpha)) * np.clip(s, 0.0, None)
    return coo_matrix((S, (i, j)), shape=W.shape).tocsr()


def build_rois(entropy, centroids, tx_xy, buffer, delta=1.0,
               entropy_quantile=0.9, budget=None):
    """Entropy-ranked ROIs; overlaps retained, ranked by rho_i (Fig. 1d).
    centroids: (N,2) as (y,x) [regionprops order]; tx_xy: (T,2) as (x,y)."""
    entropy = np.asarray(entropy, float)
    thr = float(np.quantile(entropy, entropy_quantile))
    cand = np.where(entropy >= thr)[0]
    b = float(buffer)
    txx, txy = tx_xy[:, 0], tx_xy[:, 1]
    rois = []
    for i in cand:
        y, x = centroids[i, 0], centroids[i, 1]
        x0, x1, y0, y1 = x - b, x + b, y - b, y + b
        C_i = int(np.sum((txx >= x0) & (txx <= x1) & (txy >= y0) & (txy <= y1)))
        rois.append(ROI(int(i), (x0, x1, y0, y1), float(entropy[i]), C_i,
                        float(entropy[i] / max(C_i, delta))))
    rois.sort(key=lambda r: r.utility_cost, reverse=True)
    return rois[:int(budget)] if budget else rois


def characterize_roi(roi, expr, centroids, labels_state, genes, top_n=5):
    """Summarize one ROI in place (no merging): member cells, top transcripts,
    dominant state and a state-mixing index (0 = pure, ->1 = multi-cell-type)."""
    x0, x1, y0, y1 = roi.bounds
    cy, cx = centroids[:, 0], centroids[:, 1]
    members = np.where((cx >= x0) & (cx <= x1) & (cy >= y0) & (cy <= y1))[0]
    roi.member_cells = members
    if members.size == 0:
        return {"center_cell": roi.center_cell, "entropy": roi.entropy,
                "n_transcripts": roi.n_transcripts, "utility_cost": roi.utility_cost,
                "n_members": 0, "top_transcripts": [], "dominant_state": None,
                "state_mixing": float("nan")}
    prof = np.asarray(expr, float)[members].mean(0)
    top = [str(genes[t]) for t in np.argsort(-prof)[:top_n]]
    st, cnt = np.unique(labels_state[members], return_counts=True)
    return {"center_cell": roi.center_cell, "entropy": roi.entropy,
            "n_transcripts": roi.n_transcripts, "utility_cost": roi.utility_cost,
            "n_members": int(members.size), "top_transcripts": top,
            "dominant_state": int(st[np.argmax(cnt)]),
            "state_mixing": float(1.0 - cnt.max() / cnt.sum())}


def run_discovery(adata, cluster_key, tx_xy, genes, out_dir,
                  coord_keys=("x_centroid_px", "y_centroid_px"),
                  buffer_px=50.0, entropy_quantile=0.9, interaction_alpha=1.0,
                  k=15, budget=None):
    """Diagnostic driver. Reads the coupled cell-state posterior q_i / entropy H_i
    from `adata` (written by couple_pj_to_qi), computes the Fig. 1d metrics
    (H_i, M_i, S_ij, entropy-ranked ROIs), and writes them to `out_dir`. Returns
    the summary dict. Read-only w.r.t. cell identities. Non-fatal: logs and
    returns None on any failure."""
    import pandas as pd
    try:
        os.makedirs(out_dir, exist_ok=True)
        n_obs = int(adata.n_obs)

        # entropy H_i: prefer the value written by the coupling step; else recompute
        # from the stored posterior q_coupled.
        if "posterior_entropy_coupled" in adata.obs:
            H = adata.obs["posterior_entropy_coupled"].to_numpy(float)
        elif "q_coupled" in adata.obsm:
            H = cell_state_entropy(adata.obsm["q_coupled"])
        else:
            print("[DISCOVERY] no coupled posterior/entropy found "
                  "(run couple_pj_to_qi first); skipping.", flush=True)
            return None

        # alignment guard: every per-cell array must match adata's cell count,
        # otherwise a subsampled marker fit would misalign ROIs to the wrong cells.
        if H.shape[0] != n_obs:
            print(f"[DISCOVERY] entropy length {H.shape[0]} != n_obs {n_obs}; "
                  "the coupled posterior is not aligned to adata (was the marker "
                  "fit subsampled?). Skipping to avoid misassignment.", flush=True)
            return None

        centroids = adata.obs[list(coord_keys)].to_numpy(float)[:, ::-1]  # -> (y,x)
        X = adata.X
        expr = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        labels = pd.Categorical(adata.obs[cluster_key]).codes

        # boundary score M_i needs the full posterior q; None if only entropy exists
        q = adata.obsm["q_coupled"] if "q_coupled" in adata.obsm else None
        M = boundary_score(q) if q is not None else None

        # spatial graph + interaction niches S_ij
        W = build_spatial_weights(centroids, k=int(k))
        S = interaction_score(W, centroids, expr, alpha=float(interaction_alpha))

        # entropy-ranked ROIs (overlaps retained) + in-place characterization
        rois = build_rois(H, centroids, tx_xy, buffer_px,
                          entropy_quantile=entropy_quantile, budget=budget)
        summ = [characterize_roi(r, expr, centroids, labels, genes) for r in rois]

        # ---- write per-cell diagnostics + ROI/interaction tables ----
        pd.DataFrame(summ).to_csv(os.path.join(out_dir, "discovery_rois.csv"),
                                  index=False)
        adata.obs["discovery_entropy"] = H
        adata.obs["boundary_score"] = M if M is not None else np.nan
        cols = ["discovery_entropy"] + (["boundary_score"] if M is not None else [])
        adata.obs[cols].to_csv(os.path.join(out_dir, "per_cell_uncertainty.csv"))

        Sc = S.tocoo()
        eidx = np.argsort(-Sc.data)[:min(1000, Sc.data.size)]
        pd.DataFrame({"cell_i": Sc.row[eidx], "cell_j": Sc.col[eidx],
                      "S_ij": Sc.data[eidx]}).to_csv(
            os.path.join(out_dir, "interaction_niches.csv"), index=False)

        with open(os.path.join(out_dir, "discovery_summary.json"), "w") as f:
            json.dump({"n_cells": n_obs, "n_rois": len(rois),
                       "entropy_quantile": entropy_quantile, "buffer_px": buffer_px,
                       "interaction_alpha": interaction_alpha,
                       "boundary_score_available": bool(M is not None),
                       "note": ("diagnostic module (Fig. 1d); read-only, not used for "
                                "the reported results; closed-loop refinement is a "
                                "planned extension, not implemented")}, f, indent=2)
        print(f"[DISCOVERY] wrote {len(rois)} ROIs + interaction/boundary to {out_dir}",
              flush=True)
        return {"rois": summ, "interaction": S, "boundary_score": M, "entropy": H}
    except Exception as e:
        print(f"[DISCOVERY] skipped (non-fatal): {type(e).__name__}: {e}", flush=True)
        return None
