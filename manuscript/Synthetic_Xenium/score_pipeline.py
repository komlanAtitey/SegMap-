#!/usr/bin/env python3
"""
score_pipeline.py
=================

Ground-truth scorer for spatial-transcriptomics pipelines, designed to run
against the output of ``synthetic_tissue.py``. It quantifies four tasks:

  1. SEGMENTATION        — do predicted cells match true cells?
  2. TRANSCRIPT ASSIGNMENT — are molecules assigned to the right cell?
  3. MARKER DETECTION    — are the true cell-type markers recoverable from the
                           method's cells?
  4. DENOISING           — is background (noise) separated from true signal?

Method-agnostic interface
--------------------------
The thing being scored is supplied as a predicted transcript -> cell assignment,
one row per ground-truth transcript, in the SAME ORDER as the ground-truth
``transcripts.csv``:

    a CSV/Parquet with a column `pred_cell_id`  (-1 = unassigned / predicted noise)
    optionally a column `noise_prob`            (higher = more likely background)

Optionally, a predicted integer label mask (.npy, same raster grid as the
ground-truth `label_mask_cells.npy`) enables pixel-level segmentation metrics.
If no mask is given, segmentation is scored at the object level from the
molecule-to-cell sets (always available, method-independent).

Any method (SegMap, Baysor, Cellpose+assignment, …) is scored identically by
emitting this one assignment vector, so results are directly comparable.

Outputs
-------
metrics.json (all scalars), marker_detection.csv (per cell type),
segmentation_matches.csv, and scorecard.png.

Dependencies: numpy, pandas, scipy, scikit-learn, matplotlib.
"""
from __future__ import annotations

import argparse
import json
import os
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import ranksums
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             roc_auc_score, average_precision_score)


# ============================================================
# Loading
# ============================================================

def load_ground_truth(gt_dir: str) -> dict:
    tx = pd.read_csv(os.path.join(gt_dir, "transcripts.csv"))
    panel = pd.read_csv(os.path.join(gt_dir, "gene_panel.csv"))
    cells = pd.read_csv(os.path.join(gt_dir, "cells.csv"))
    with open(os.path.join(gt_dir, "ground_truth_summary.json")) as f:
        summary = json.load(f)
    n_genes = len(panel)
    # true marker sets per type
    true_markers = {int(t): set(panel.loc[panel.marker_of_type == t, "gene_index"].tolist())
                    for t in sorted(panel.loc[panel.marker_of_type >= 0, "marker_of_type"].unique())}
    gt_mask_path = os.path.join(gt_dir, "label_mask_cells.npy")
    gt_mask = np.load(gt_mask_path) if os.path.isfile(gt_mask_path) else None
    return dict(tx=tx, panel=panel, cells=cells, summary=summary,
                n_genes=int(n_genes), true_markers=true_markers, gt_mask=gt_mask)


def load_prediction(pred_path: str, n_tx: int,
                    noise_col: str = "noise_prob") -> Tuple[np.ndarray, Optional[np.ndarray]]:
    if pred_path.endswith(".parquet"):
        pred = pd.read_parquet(pred_path)
    else:
        pred = pd.read_csv(pred_path)
    if "pred_cell_id" not in pred.columns:
        raise ValueError("Prediction file must contain a 'pred_cell_id' column.")
    if "index" in pred.columns and len(pred) == n_tx:
        pred = pred.sort_values("index")
    if len(pred) != n_tx:
        raise ValueError(f"Prediction has {len(pred)} rows but ground truth has "
                         f"{n_tx} transcripts; they must align row-for-row.")
    pred_cell = pred["pred_cell_id"].to_numpy()
    noise_prob = pred[noise_col].to_numpy() if noise_col in pred.columns else None
    return pred_cell, noise_prob


# ============================================================
# Cell matching (molecule-set IoU + Hungarian)
# ============================================================

def _contingency(true_lab: np.ndarray, pred_lab: np.ndarray):
    """Dense intersection matrix over molecules where both labels are >= 0."""
    t_ids = np.unique(true_lab[true_lab >= 0])
    p_ids = np.unique(pred_lab[pred_lab >= 0])
    t_pos = {t: i for i, t in enumerate(t_ids)}
    p_pos = {p: j for j, p in enumerate(p_ids)}
    inter = np.zeros((len(t_ids), len(p_ids)), dtype=np.int64)
    m = (true_lab >= 0) & (pred_lab >= 0)
    ti = np.array([t_pos[t] for t in true_lab[m]])
    pj = np.array([p_pos[p] for p in pred_lab[m]])
    if len(ti):
        np.add.at(inter, (ti, pj), 1)
    t_size = np.array([(true_lab == t).sum() for t in t_ids])
    p_size = np.array([(pred_lab == p).sum() for p in p_ids])
    return t_ids, p_ids, t_size, p_size, inter


def match_cells(true_lab: np.ndarray, pred_lab: np.ndarray):
    t_ids, p_ids, t_size, p_size, inter = _contingency(true_lab, pred_lab)
    if len(t_ids) == 0 or len(p_ids) == 0:
        return dict(t_ids=t_ids, p_ids=p_ids, t_size=t_size, p_size=p_size,
                    inter=inter, iou=np.zeros((len(t_ids), len(p_ids))),
                    row=np.array([], int), col=np.array([], int))
    union = t_size[:, None] + p_size[None, :] - inter
    iou = np.where(union > 0, inter / union, 0.0)
    row, col = linear_sum_assignment(-iou)           # maximize total IoU
    keep = iou[row, col] > 0
    return dict(t_ids=t_ids, p_ids=p_ids, t_size=t_size, p_size=p_size,
                inter=inter, iou=iou, row=row[keep], col=col[keep])


# ============================================================
# 1. Segmentation
# ============================================================

def score_segmentation(match: dict, thresholds=(0.3, 0.5, 0.7)) -> dict:
    iou = match["iou"]; row, col = match["row"], match["col"]
    n_true, n_pred = len(match["t_ids"]), len(match["p_ids"])
    matched_iou = iou[row, col] if len(row) else np.array([])
    out = {"n_true_cells": int(n_true), "n_pred_cells": int(n_pred),
           "count_error": int(n_pred - n_true),
           "mean_matched_iou": float(np.mean(matched_iou)) if matched_iou.size else 0.0,
           "median_matched_iou": float(np.median(matched_iou)) if matched_iou.size else 0.0}
    for tau in thresholds:
        tp = int(np.sum(matched_iou >= tau))
        fp = n_pred - tp
        fn = n_true - tp
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn) if (tp + fn) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        out[f"detF1@{tau}"] = round(f1, 4)
        out[f"precision@{tau}"] = round(prec, 4)
        out[f"recall@{tau}"] = round(rec, 4)
    # split / merge: coverage >= 20% of a cell's molecules
    inter = match["inter"]; tsz = match["t_size"]; psz = match["p_size"]
    if inter.size:
        cov_true = inter / np.maximum(tsz[:, None], 1)     # fraction of true cell in pred
        cov_pred = inter / np.maximum(psz[None, :], 1)     # fraction of pred cell from true
        out["over_segmentation_rate"] = float(np.mean(np.sum(cov_true >= 0.2, axis=1) > 1))
        out["under_segmentation_rate"] = float(np.mean(np.sum(cov_pred >= 0.2, axis=0) > 1))
    else:
        out["over_segmentation_rate"] = out["under_segmentation_rate"] = 0.0
    return out


def score_segmentation_pixels(gt_mask: np.ndarray, pred_mask: np.ndarray,
                              thresholds=(0.5,)) -> dict:
    if gt_mask.shape != pred_mask.shape:
        raise ValueError(f"mask shapes differ: {gt_mask.shape} vs {pred_mask.shape}")
    tl = gt_mask.ravel(); pl = pred_mask.ravel()
    match = match_cells(tl.astype(np.int64) - 1, pl.astype(np.int64) - 1)  # 0=bg -> -1
    seg = score_segmentation(match, thresholds=thresholds)
    return {f"pixel_{k}": v for k, v in seg.items()}


# ============================================================
# 2. Transcript assignment
# ============================================================

def build_cell_gene(labels: np.ndarray, gene_index: np.ndarray, n_genes: int):
    """Cell x gene count matrix from a molecule->cell labelling (labels>=0)."""
    ids = np.unique(labels[labels >= 0])
    pos = {c: i for i, c in enumerate(ids)}
    M = np.zeros((len(ids), n_genes), dtype=np.float64)
    m = labels >= 0
    ri = np.array([pos[c] for c in labels[m]])
    gi = gene_index[m].astype(int)
    if len(ri):
        np.add.at(M, (ri, gi), 1)
    return ids, M


def score_assignment(true_lab: np.ndarray, pred_lab: np.ndarray, match: dict,
                     gene_index: np.ndarray, n_genes: int) -> dict:
    intra = true_lab >= 0                      # molecules with a true home
    # matched accuracy: true cell t -> its matched pred cell
    t_ids, p_ids = match["t_ids"], match["p_ids"]
    t2p = {int(t_ids[r]): int(p_ids[c]) for r, c in zip(match["row"], match["col"])}
    expected_pred = np.array([t2p.get(int(t), -999) for t in true_lab[intra]])
    acc = float(np.mean(pred_lab[intra] == expected_pred))

    # ARI/NMI over true intracellular molecules. Unassigned (pred<0) molecules
    # are treated as SINGLETONS, not one shared cluster: a dropped molecule is
    # "not grouped", not "mis-merged with every other dropped molecule". Lumping
    # all -1 together creates a huge false-merge that collapses ARI even when
    # per-molecule accuracy and purity are high. Signal retention (below) and
    # completeness separately capture how much was dropped.
    pl_fair = pred_lab[intra].astype(np.int64).copy()
    un = pl_fair < 0
    if un.any():
        pl_fair[un] = -(np.arange(1, int(un.sum()) + 1))   # unique singleton ids
    ari = adjusted_rand_score(true_lab[intra], pl_fair)
    nmi = normalized_mutual_info_score(true_lab[intra], pl_fair)
    # transparency: the old (unassigned-as-one-cluster) value, for reference
    ari_incl = adjusted_rand_score(true_lab[intra], pred_lab[intra])

    # purity / completeness from contingency
    inter = match["inter"]; tsz = match["t_size"]; psz = match["p_size"]
    purity = float(np.mean(inter.max(axis=0) / np.maximum(psz, 1))) if inter.size else 0.0
    if len(match["row"]):
        comp = inter[match["row"], match["col"]] / np.maximum(tsz[match["row"]], 1)
        completeness = float(np.mean(comp))
    else:
        completeness = 0.0

    # expression fidelity on matched cells
    t_ids2, Mt = build_cell_gene(true_lab, gene_index, n_genes)
    p_ids2, Mp = build_cell_gene(pred_lab, gene_index, n_genes)
    tpos = {int(c): i for i, c in enumerate(t_ids2)}
    ppos = {int(c): i for i, c in enumerate(p_ids2)}
    pear, spear = [], []
    for r, c in zip(match["row"], match["col"]):
        tt, pp = int(t_ids[r]), int(p_ids[c])
        if tt in tpos and pp in ppos:
            a = np.log1p(Mt[tpos[tt]]); b = np.log1p(Mp[ppos[pp]])
            if a.std() > 0 and b.std() > 0:
                pear.append(np.corrcoef(a, b)[0, 1])
                ra = pd.Series(Mt[tpos[tt]]).rank().to_numpy()
                rb = pd.Series(Mp[ppos[pp]]).rank().to_numpy()
                spear.append(np.corrcoef(ra, rb)[0, 1])
    return {"assignment_accuracy": round(acc, 4),
            "ARI": round(float(ari), 4),
            "ARI_incl_unassigned": round(float(ari_incl), 4),
            "NMI": round(float(nmi), 4),
            "mean_cell_purity": round(purity, 4),
            "mean_cell_completeness": round(completeness, 4),
            "expr_pearson_log_mean": round(float(np.mean(pear)), 4) if pear else 0.0,
            "expr_spearman_mean": round(float(np.mean(spear)), 4) if spear else 0.0}


# ============================================================
# 3. Marker detection
# ============================================================

def _bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def _markers_from_cells(M: np.ndarray, cell_types: np.ndarray, n_genes: int,
                        padj=0.05, logfc=0.25, min_pct=0.25):
    """One-vs-rest Wilcoxon marker calling across cells. Returns {type: (R_set, scores)}."""
    lib = M.sum(1, keepdims=True); lib[lib == 0] = 1
    norm = M / lib * 1e4
    log = np.log1p(norm)
    out = {}
    types = [t for t in np.unique(cell_types) if t >= 0]
    for t in types:
        grp = cell_types == t; rest = ~grp
        if grp.sum() < 2 or rest.sum() < 2:
            out[int(t)] = (set(), np.zeros(n_genes)); continue
        pvals = np.ones(n_genes); scores = np.zeros(n_genes)
        for g in range(n_genes):
            a, b = log[grp, g], log[rest, g]
            if a.any() or b.any():
                s, p = ranksums(a, b); scores[g] = s; pvals[g] = p
        q = _bh(pvals)
        mg = norm[grp, :].mean(0); mr = norm[rest, :].mean(0)
        lfc = np.log2((mg + 1e-9) / (mr + 1e-9))
        pct = (M[grp, :] > 0).mean(0)
        R = set(np.where((q < padj) & (lfc >= logfc) & (pct >= min_pct))[0].tolist())
        out[int(t)] = (R, scores)
    return out


def score_markers(true_lab, pred_lab, true_cell_type, gene_index, n_genes,
                  true_markers: Dict[int, set]) -> dict:
    """Marker recoverability: build predicted cells, label each by the majority
    TRUE type of its molecules, then compare recovered markers to the true sets.
    Reported alongside the ceiling computed on the true cells."""
    def typed_matrix(labels):
        ids, M = build_cell_gene(labels, gene_index, n_genes)
        # majority true type per cell
        ctype = np.full(len(ids), -1)
        pos = {int(c): i for i, c in enumerate(ids)}
        df = pd.DataFrame({"lab": labels, "typ": true_cell_type})
        df = df[(df.lab >= 0) & (df.typ >= 0)]
        if len(df):
            maj = df.groupby("lab")["typ"].agg(lambda s: s.value_counts().index[0])
            for c, t in maj.items():
                if int(c) in pos:
                    ctype[pos[int(c)]] = int(t)
        return M, ctype

    def prf(recovered):
        rows = []
        for t, M_t in true_markers.items():
            R = recovered.get(t, (set(), None))[0]
            tp = len(R & M_t)
            prec = tp / len(R) if R else 0.0
            rec = tp / len(M_t) if M_t else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            # average precision from DE scores
            sc = recovered.get(t, (None, np.zeros(n_genes)))[1]
            y = np.zeros(n_genes); y[list(M_t)] = 1
            ap = average_precision_score(y, np.abs(sc)) if sc is not None and y.sum() else 0.0
            rows.append(dict(cell_type=t, n_true_markers=len(M_t), n_recovered=len(R),
                             precision=round(prec, 4), recall=round(rec, 4),
                             f1=round(f1, 4), average_precision=round(float(ap), 4)))
        return pd.DataFrame(rows)

    Mp, tp_type = typed_matrix(pred_lab)
    Mt, tt_type = typed_matrix(true_lab)
    rec_pred = _markers_from_cells(Mp, tp_type, n_genes)
    rec_true = _markers_from_cells(Mt, tt_type, n_genes)
    df_pred = prf(rec_pred); df_true = prf(rec_true)
    summary = {"marker_f1_mean": round(float(df_pred.f1.mean()), 4),
               "marker_precision_mean": round(float(df_pred.precision.mean()), 4),
               "marker_recall_mean": round(float(df_pred.recall.mean()), 4),
               "marker_AP_mean": round(float(df_pred.average_precision.mean()), 4),
               "marker_f1_ceiling_true_cells": round(float(df_true.f1.mean()), 4)}
    return summary, df_pred


# ============================================================
# 4. Denoising
# ============================================================

def score_denoising(true_lab, pred_lab, noise_prob=None) -> dict:
    is_noise = (true_lab < 0)                       # positive class = background
    pred_noise = (pred_lab < 0)
    tp = int(np.sum(is_noise & pred_noise))
    fp = int(np.sum(~is_noise & pred_noise))        # real signal wrongly dropped
    fn = int(np.sum(is_noise & ~pred_noise))        # noise wrongly kept
    tn = int(np.sum(~is_noise & ~pred_noise))
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    out = {"noise_detection_precision": round(prec, 4),
           "noise_detection_recall": round(rec, 4),
           "noise_detection_f1": round(f1, 4),
           "signal_retention": round(tn / (tn + fp) if (tn + fp) else 0.0, 4),
           "signal_purity": round(tn / (tn + fn) if (tn + fn) else 0.0, 4)}
    if noise_prob is not None and is_noise.any() and (~is_noise).any():
        out["noise_roc_auc"] = round(float(roc_auc_score(is_noise, noise_prob)), 4)
        out["noise_pr_auc"] = round(float(average_precision_score(is_noise, noise_prob)), 4)
    return out


# ============================================================
# Orchestration
# ============================================================

def score_all(gt: dict, pred_cell: np.ndarray, noise_prob=None,
              pred_mask: Optional[np.ndarray] = None) -> Tuple[dict, pd.DataFrame]:
    tx = gt["tx"]
    true_lab = tx["cell_id"].to_numpy().astype(np.int64)
    gene_index = tx["gene_index"].to_numpy().astype(np.int64)
    true_type = tx["cell_type"].to_numpy().astype(np.int64)
    pred_lab = np.asarray(pred_cell).astype(np.int64)
    n_genes = gt["n_genes"]

    match = match_cells(true_lab, pred_lab)
    metrics = {"segmentation": score_segmentation(match),
               "assignment": score_assignment(true_lab, pred_lab, match,
                                               gene_index, n_genes),
               "denoising": score_denoising(true_lab, pred_lab, noise_prob)}
    marker_summary, marker_df = score_markers(true_lab, pred_lab, true_type,
                                              gene_index, n_genes, gt["true_markers"])
    metrics["markers"] = marker_summary
    if pred_mask is not None and gt["gt_mask"] is not None:
        metrics["segmentation"].update(
            score_segmentation_pixels(gt["gt_mask"], pred_mask))
    return metrics, marker_df


def scorecard_figure(metrics: dict, path: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [
        ("Seg detF1@0.5", metrics["segmentation"].get("detF1@0.5", 0)),
        ("Seg mean IoU", metrics["segmentation"].get("mean_matched_iou", 0)),
        ("Assign acc", metrics["assignment"]["assignment_accuracy"]),
        ("Assign ARI", metrics["assignment"]["ARI"]),
        ("Cell purity", metrics["assignment"]["mean_cell_purity"]),
        ("Expr r (log)", metrics["assignment"]["expr_pearson_log_mean"]),
        ("Marker F1", metrics["markers"]["marker_f1_mean"]),
        ("Marker AP", metrics["markers"]["marker_AP_mean"]),
        ("Denoise F1", metrics["denoising"]["noise_detection_f1"]),
        ("Signal keep", metrics["denoising"]["signal_retention"]),
    ]
    names = [a for a, _ in labels]; vals = [max(0, float(b)) for _, b in labels]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    bars = ax.bar(range(len(vals)), vals, color="#1f77b4")
    ax.set_xticks(range(len(vals))); ax.set_xticklabels(names, rotation=35, ha="right")
    ax.set_ylim(0, 1.05); ax.set_ylabel("score"); ax.axhline(1.0, color="0.7", lw=0.8, ls="--")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                ha="center", va="bottom", fontsize=8)
    ax.set_title("Benchmark scorecard (higher = better; 1.0 = perfect)", loc="left")
    plt.tight_layout(); fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


# ============================================================
# CLI
# ============================================================

def build_parser():
    p = argparse.ArgumentParser(description="Score a pipeline against synthetic "
                                            "ground truth (segmentation, assignment, "
                                            "markers, denoising).")
    p.add_argument("--gt-dir", required=True, help="synthetic_tissue.py output dir.")
    p.add_argument("--pred", required=True,
                   help="CSV/Parquet with 'pred_cell_id' (row-aligned to transcripts.csv).")
    p.add_argument("--pred-mask", default=None, help="Optional predicted label .npy mask.")
    p.add_argument("--noise-col", default="noise_prob", help="Optional noise-score column.")
    p.add_argument("--out-dir", default="score_out")
    p.add_argument("--no-figure", action="store_true")
    return p


def main():
    a = build_parser().parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    gt = load_ground_truth(a.gt_dir)
    pred_cell, noise_prob = load_prediction(a.pred, len(gt["tx"]), a.noise_col)
    pred_mask = np.load(a.pred_mask) if a.pred_mask else None

    metrics, marker_df = score_all(gt, pred_cell, noise_prob, pred_mask)

    with open(os.path.join(a.out_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2, default=float)
    marker_df.to_csv(os.path.join(a.out_dir, "marker_detection.csv"), index=False)
    if not a.no_figure:
        scorecard_figure(metrics, os.path.join(a.out_dir, "scorecard.png"))

    print("\n================ BENCHMARK METRICS ================")
    for task, d in metrics.items():
        print(f"\n[{task}]")
        for k, v in d.items():
            print(f"  {k}: {v}")
    print(f"\nWritten to: {os.path.abspath(a.out_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
