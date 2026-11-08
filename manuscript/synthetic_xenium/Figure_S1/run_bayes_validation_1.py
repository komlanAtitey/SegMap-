#!/usr/bin/env python3
# =============================================================================
# run_bayes_validation.py  --  COMPLETE, SELF-CONTAINED validation script.
#
# Produces the four validation plots for the SegMap Bayesian marker model AND
# writes the underlying data table behind each plot, so you can inspect/collect
# every number that goes into a figure.
#
#   Figures (in ./bayes_validation_out/figures/):
#     roc.png                  ROC of marker detection            (needs truth -> simulation)
#     fsv_recovery.png         estimated vs true FSV              (needs truth -> simulation)
#     null_calibration.png     posteriors under permuted coords   (calibration)
#     posterior_separation.png posterior by true spatial label
#
#   Data collected (in ./bayes_validation_out/data/):
#     per_gene_fit.csv         gene, posterior, fsv, logBF, is_spatial, true_fsv
#     roc_curve.csv            fpr, tpr  (+ AUC in summary)
#     fsv_recovery.csv         gene, true_fsv, est_fsv
#     null_posteriors.csv      perm, gene, posterior
#     summary.csv              AUC, FSV_RMSE, null_pi_hat, n_selected, precision, recall
#
# RUN:   python run_bayes_validation.py
# DEPS:  numpy, scipy, scikit-learn, matplotlib   (all standard)
# =============================================================================
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, diags, identity
from scipy.sparse.linalg import eigsh
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------- CONFIG ------------------------------------- #
N, G, FRAC_SPATIAL, LENGTHSCALE = 1500, 300, 0.30, 0.10   # simulation size / structure
FSV_RANGE = (0.30, 0.80)                                  # true FSV of spatial genes
K_NN, M_BASIS = 15, 100                                   # graph k, spatial basis size
ALPHA = 0.05                                              # Bayesian FDR level
N_PERM = 10                                               # null-calibration permutations
SEED = 1
OUTDIR = "bayes_validation_out"


# ===================== Bayesian marker model (inline) ====================== #
def build_spatial_basis(coords, k=15, m=100, normalized=True):
    coords = np.asarray(coords, float); Nn = coords.shape[0]
    m = int(min(m, max(1, Nn - 2))); k = int(min(k, max(1, Nn - 1)))
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(coords)
    idx = idx[:, 1:]
    rows = np.repeat(np.arange(Nn), k); cols = idx.reshape(-1)
    W = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(Nn, Nn))
    W = ((W + W.T) > 0).astype(float)
    deg = np.asarray(W.sum(1)).ravel()
    if normalized:
        dinv = 1.0 / np.sqrt(np.maximum(deg, 1e-8))
        L = identity(Nn) - diags(dinv) @ W @ diags(dinv)
    else:
        L = diags(deg) - W
    vals, vecs = eigsh(L.tocsc(), k=min(m + 1, Nn - 1), sigma=-1e-6, which="LM")
    o = np.argsort(vals); vals, vecs = vals[o], vecs[:, o]
    keep = np.where(vals > 1e-8)[0][:m]
    return {"U": vecs[:, keep], "gamma": vals[keep], "N": Nn, "m": len(keep)}


def _nll(par, yt2, gamma, ss_perp, Nn, m):
    s2, e2 = np.exp(par)
    d = np.maximum(s2 / gamma + e2, 1e-10)
    return 0.5 * (np.sum(yt2 / d + np.log(d)) + ss_perp / e2
                  + (Nn - m) * np.log(e2) + Nn * np.log(2 * np.pi))


def _fit_gene(y, b):
    U, gamma, Nn, m = b["U"], b["gamma"], b["N"], b["m"]
    yc = y - y.mean(); yt = U.T @ yc; yt2 = yt * yt
    ss_perp = max(float(yc @ yc) - float(yt2.sum()), 0.0)
    v = max(float(yc @ yc) / Nn, 1e-8)
    opt = minimize(_nll, np.log([v / 2, v / 2]),
                   args=(yt2, gamma, ss_perp, Nn, m), method="L-BFGS-B")
    s2, e2 = np.exp(opt.x)
    return -opt.fun, -0.5 * Nn * (np.log(2 * np.pi * v) + 1.0), s2, e2, s2 / (s2 + e2)


def fit_all_genes(expr, b, dpar=1, pi0=0.5, eb_iter=100, eb_tol=1e-4):
    expr = np.asarray(expr, float); Gn = expr.shape[1]; Nn = b["N"]
    logL1 = np.empty(Gn); logL0 = np.empty(Gn); fsv = np.empty(Gn)
    for j in range(Gn):
        logL1[j], logL0[j], _, _, fsv[j] = _fit_gene(expr[:, j], b)
    logBF = logL1 - logL0 - 0.5 * dpar * np.log(Nn)
    sig = lambda x: 1.0 / (1.0 + np.exp(-x)); logit = lambda p: np.log(p / (1 - p))
    pi = pi0
    for _ in range(eb_iter):
        new = min(max(float(sig(logBF + logit(pi)).mean()), 1e-4), 1 - 1e-4)
        if abs(new - pi) < eb_tol:
            pi = new; break
        pi = new
    return {"logBF": logBF, "posterior": sig(logBF + logit(pi)), "fsv": fsv, "pi_hat": pi}


def select_markers_bfdr(posterior, alpha=0.05):
    o = np.argsort(-posterior)
    fdr = np.cumsum(1.0 - posterior[o]) / (np.arange(o.size) + 1)
    keep = np.where(fdr <= alpha)[0]
    return o[: keep.max() + 1] if keep.size else np.array([], dtype=int)


# =============================== simulation ================================ #
def simulate(N, G, frac, lengthscale, fsv_range, seed):
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2))
    D2 = ((coords[:, None, :] - coords[None, :, :]) ** 2).sum(-1)
    Kc = np.exp(-D2 / (2 * lengthscale ** 2)) + 1e-6 * np.eye(N)
    L = np.linalg.cholesky(Kc)
    n_sp = round(frac * G); is_sp = np.zeros(G, bool); is_sp[:n_sp] = True
    fsv_true = np.zeros(G); expr = np.empty((N, G))
    for j in range(G):
        if is_sp[j]:
            f = rng.uniform(*fsv_range); fsv_true[j] = f
            expr[:, j] = np.sqrt(f) * (L @ rng.standard_normal(N)) + \
                rng.normal(0, np.sqrt(1 - f), N)
        else:
            expr[:, j] = rng.standard_normal(N)
    genes = [f"g{j+1}" for j in range(G)]
    return coords, expr, is_sp, fsv_true, genes


def roc_points(score, label):
    o = np.argsort(-score); y = label[o]
    P = max(y.sum(), 1); Nn = max((~label).sum(), 1)
    return (np.r_[0, np.cumsum(~y) / Nn], np.r_[0, np.cumsum(y) / P],
            float(np.trapezoid(np.r_[0, np.cumsum(y) / P], np.r_[0, np.cumsum(~y) / Nn])))


# ================================= main ==================================== #
def main():
    fig_dir = os.path.join(OUTDIR, "figures"); data_dir = os.path.join(OUTDIR, "data")
    os.makedirs(fig_dir, exist_ok=True); os.makedirs(data_dir, exist_ok=True)

    print("[1/4] simulating data with known ground truth ...")
    coords, expr, is_sp, fsv_true, genes = simulate(
        N, G, FRAC_SPATIAL, LENGTHSCALE, FSV_RANGE, SEED)

    print("[2/4] fitting Bayesian marker model ...")
    basis = build_spatial_basis(coords, k=K_NN, m=M_BASIS)
    fit = fit_all_genes(expr, basis)
    sel = select_markers_bfdr(fit["posterior"], ALPHA)

    # ---- collect per-gene data ----
    per_gene = pd.DataFrame({
        "gene": genes, "posterior": fit["posterior"], "fsv": fit["fsv"],
        "logBF": fit["logBF"], "is_spatial": is_sp.astype(int), "true_fsv": fsv_true,
    }).sort_values("posterior", ascending=False)
    per_gene.to_csv(os.path.join(data_dir, "per_gene_fit.csv"), index=False)

    # ---- ROC data + figure ----
    fpr, tpr, auc = roc_points(fit["posterior"], is_sp)
    pd.DataFrame({"fpr": fpr, "tpr": tpr}).to_csv(
        os.path.join(data_dir, "roc_curve.csv"), index=False)
    plt.figure(figsize=(4, 4))
    plt.plot(fpr, tpr, color="#2c7fb8", lw=2); plt.plot([0, 1], [0, 1], "--", color="grey")
    plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
    plt.title(f"Marker detection ROC (AUC = {auc:.3f})")
    plt.gca().set_aspect("equal"); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "roc.png"), dpi=150); plt.close()

    # ---- FSV recovery data + figure ----
    t, e = fsv_true[is_sp], fit["fsv"][is_sp]
    fsv_rmse = float(np.sqrt(np.mean((e - t) ** 2)))
    pd.DataFrame({"gene": np.array(genes)[is_sp], "true_fsv": t, "est_fsv": e}).to_csv(
        os.path.join(data_dir, "fsv_recovery.csv"), index=False)
    plt.figure(figsize=(4, 4))
    plt.scatter(t, e, s=12, alpha=0.6, color="#756bb1"); plt.plot([0, 1], [0, 1], "--", color="grey")
    plt.xlabel("True FSV"); plt.ylabel("Estimated FSV")
    plt.title(f"FSV recovery (RMSE = {fsv_rmse:.3f})")
    plt.gca().set_aspect("equal"); plt.xlim(0, 1); plt.ylim(0, 1); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "fsv_recovery.png"), dpi=150); plt.close()

    # ---- posterior separation figure ----
    plt.figure(figsize=(4, 4))
    plt.violinplot([fit["posterior"][~is_sp], fit["posterior"][is_sp]], showmedians=True)
    plt.xticks([1, 2], ["null", "spatial"]); plt.ylabel("Posterior P(X_j = 1)")
    plt.title("Posterior by true label"); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "posterior_separation.png"), dpi=150); plt.close()

    # ---- null calibration: permute coordinates, collect posteriors ----
    print(f"[3/4] null calibration ({N_PERM} coordinate permutations) ...")
    rng = np.random.default_rng(SEED + 100); rows = []; pis = []
    for b in range(N_PERM):
        perm = coords[rng.permutation(N)]
        fb = fit_all_genes(expr, build_spatial_basis(perm, k=K_NN, m=M_BASIS))
        pis.append(fb["pi_hat"])
        rows.append(pd.DataFrame({"perm": b, "gene": genes, "posterior": fb["posterior"]}))
    null_df = pd.concat(rows, ignore_index=True)
    null_df.to_csv(os.path.join(data_dir, "null_posteriors.csv"), index=False)
    pih = float(np.mean(pis))
    plt.figure(figsize=(4.5, 3.5))
    plt.hist(null_df["posterior"], bins=40, color="#bdbdbd", edgecolor="white")
    plt.axvline(pih, ls="--", color="#2c7fb8")
    plt.xlabel("Posterior under permuted coordinates"); plt.ylabel("Gene count")
    plt.title(f"Null calibration (pi_hat = {pih:.3f})"); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "null_calibration.png"), dpi=150); plt.close()

    # ---- summary ----
    prec = float(is_sp[sel].mean()) if sel.size else float("nan")
    rec = float(is_sp[sel].sum() / is_sp.sum())
    summary = pd.DataFrame([{
        "AUC": auc, "FSV_RMSE": fsv_rmse, "null_pi_hat": pih,
        "n_selected": int(sel.size), "precision": prec, "recall": rec,
        "pi_hat_spatial": fit["pi_hat"], "N": N, "G": G,
    }])
    summary.to_csv(os.path.join(data_dir, "summary.csv"), index=False)

    print("[4/4] done.")
    print(summary.to_string(index=False))
    print(f"\nFigures -> {fig_dir}/    Data -> {data_dir}/")


if __name__ == "__main__":
    main()
