#!/usr/bin/env python3
"""
segmap_to_pred.py  —  OUTPUT adapter (SegMap cell mask -> scorer-ready prediction)

SegMap writes its final expanded cell mask to
  <outdir>/data/full_mask_cellpose_fullres_tiled.tif
but never saves a per-transcript assignment. This reconstructs the prediction
the way SegMap assigns transcripts internally: look up each synthetic
transcript's pixel in that mask. The mask label becomes `pred_cell_id`
(0 -> -1 = unassigned / predicted noise). A distance-to-nearest-cell value is
written as `noise_prob` for threshold-free denoising AUC.

SCALE IS INFERRED AUTOMATICALLY from the mask dimensions and the ground-truth
tissue domain (scale = mask_width / domain_width), so a wrong hand-passed
--scale can no longer silently scramble the mapping. Pass --scale only to
override; a mismatch with the inferred value is reported.

Example
-------
python segmap_to_pred.py --gt-dir synthetic_out \
    --mask synth_xenium/segmap_out/data/full_mask_cellpose_fullres_tiled.tif \
    --out segmap_pred.csv
"""
from __future__ import annotations
import argparse, json, os
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import distance_transform_edt


def main():
    ap = argparse.ArgumentParser(description="Convert a SegMap cell mask to a prediction file.")
    ap.add_argument("--gt-dir", required=True, help="synthetic_tissue.py output dir.")
    ap.add_argument("--mask", required=True, help="SegMap full_mask_cellpose_fullres_tiled.tif")
    ap.add_argument("--scale", type=float, default=None,
                    help="Override pixels/unit. Default: inferred from mask size "
                         "and tissue domain (recommended).")
    ap.add_argument("--out", default="segmap_pred.csv")
    a = ap.parse_args()

    tx = pd.read_csv(os.path.join(a.gt_dir, "transcripts.csv"))
    with open(os.path.join(a.gt_dir, "ground_truth_summary.json")) as f:
        summary = json.load(f)
    domW, domH = summary["stats"]["domain"]

    mask = tifffile.imread(a.mask)
    if mask.ndim != 2:
        mask = np.squeeze(mask)
    mask = mask.astype(np.int64)
    H, W = mask.shape

    # ---- infer scale from geometry (mask pixels per tissue unit) ----
    sx, sy = W / domW, H / domH
    inferred = float(np.round((sx + sy) / 2.0, 4))
    if abs(sx - sy) / max(sx, sy) > 0.02:
        print(f"[WARN] non-uniform inferred scale: x={sx:.3f} y={sy:.3f} "
              f"(mask {W}x{H}, domain {domW:.1f}x{domH:.1f}) — image may be transposed.")
    scale = a.scale if a.scale is not None else inferred
    if a.scale is not None and abs(a.scale - inferred) / inferred > 0.02:
        print(f"[WARN] --scale {a.scale} disagrees with inferred {inferred:.3f} "
              f"from the mask; using --scale as requested. If scores look wrong, "
              f"drop --scale to use the inferred value.")
    print(f"mask {W}x{H} | domain {domW:.1f}x{domH:.1f} | scale used = {scale:.4f}"
          + ("" if a.scale is not None else " (inferred)"))

    # ---- map each transcript to its mask pixel ----
    xr = tx["x"].to_numpy() * scale
    yr = tx["y"].to_numpy() * scale
    out_of_bounds = (xr < 0) | (xr > W - 1) | (yr < 0) | (yr > H - 1)
    if out_of_bounds.mean() > 0.02:
        print(f"[WARN] {out_of_bounds.mean():.1%} of transcripts fall OUTSIDE the "
              f"mask bounds — almost certainly a scale mismatch. Re-run without "
              f"--scale to use the inferred value ({inferred:.3f}).")
    col = np.clip(np.round(xr).astype(int), 0, W - 1)
    row = np.clip(np.round(yr).astype(int), 0, H - 1)
    label = mask[row, col]
    pred = np.where(label > 0, label, -1).astype(int)

    # soft noise score: distance (px) to nearest cell pixel, monotone -> [0,1]
    dist = distance_transform_edt(mask == 0)
    noise_prob = 1.0 - np.exp(-dist[row, col] / max(scale, 1.0))

    pd.DataFrame({"pred_cell_id": pred, "noise_prob": noise_prob}).to_csv(a.out, index=False)
    n_assigned = int((pred >= 0).sum())
    print(f"transcripts={len(pred):,}  assigned={n_assigned:,} ({n_assigned/len(pred):.1%})  "
          f"unassigned/noise={len(pred)-n_assigned:,}")
    print(f"predicted cells (unique)={int(len(np.unique(pred[pred>=0])))}")
    if n_assigned / len(pred) < 0.5:
        print("[WARN] <50% of transcripts assigned. Check that --mask is the "
              "EXPANDED cell mask (full_mask_cellpose_fullres_tiled.tif), not the "
              "nucleus-only mask (full_mask_nuclei.tif).")
    print(f"wrote: {os.path.abspath(a.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

