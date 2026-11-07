#!/usr/bin/env python3
"""
synth_to_xenium.py  —  INPUT adapter (synthetic tissue -> SegMap-ready Xenium folder)

SegMap reads `--xenium_dir` and uses:
  * morphology.ome.tif  -> Cellpose nucleus segmentation (image ALONE builds the mask)
  * transcripts.parquet -> x/y/gene for assignment + QC

This renders a synthetic DAPI image from the ground-truth NUCLEUS polygons and
writes transcripts in Xenium schema, in a single pixel coordinate frame
(`--scale` pixels per synthetic length unit), so transcripts and image are
automatically co-registered (SegMap errors out if they are not).

After running SegMap on the folder this produces, use `segmap_to_pred.py` to
turn SegMap's exported cell mask into a scorer-ready prediction.

Example
-------
python synth_to_xenium.py --gt-dir synthetic_out --out-dir synth_xenium --scale 3
# -> prints the recommended --cellpose_diameter for this scale
"""
from __future__ import annotations
import argparse, json, os
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import distance_transform_edt, gaussian_filter


def _points_in_poly(pts, poly):
    x, y = pts[:, 0], pts[:, 1]
    vx, vy = poly[:, 0], poly[:, 1]
    n = len(poly); inside = np.zeros(len(pts), bool); j = n - 1
    for i in range(n):
        xi, yi, xj, yj = vx[i], vy[i], vx[j], vy[j]
        cond = ((yi > y) != (yj > y)) & (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)
        inside ^= cond; j = i
    return inside


def rasterize(polys, H, W, scale):
    """Label image (0=bg, k+1 for polys[k]) at the given pixel scale."""
    lab = np.zeros((H, W), np.int32)
    for k, poly in polys.items():
        p = np.asarray(poly) * scale
        mnx, mny = np.floor(p.min(0)).astype(int); mxx, mxy = np.ceil(p.max(0)).astype(int)
        i0, i1 = max(0, mnx), min(W, mxx); j0, j1 = max(0, mny), min(H, mxy)
        if i1 <= i0 or j1 <= j0:
            continue
        xs = np.arange(i0, i1) + 0.5; ys = np.arange(j0, j1) + 0.5
        gx, gy = np.meshgrid(xs, ys)
        inside = _points_in_poly(np.column_stack([gx.ravel(), gy.ravel()]), p).reshape(gy.shape)
        lab[j0:j1, i0:i1][inside] = k + 1
    return lab


def render_dapi(nuclei, H, W, scale, seed=0):
    """DAPI-like 16-bit image: nuclei as centre-bright blobs on a dark background."""
    rng = np.random.default_rng(seed)
    nlab = rasterize(nuclei, H, W, scale)
    binary = nlab > 0
    # centre-bright: distance to nucleus edge, normalized per image
    edt = distance_transform_edt(binary)
    inten = np.zeros((H, W), np.float32)
    if edt.max() > 0:
        inten[binary] = 0.35 + 0.65 * (edt[binary] / edt.max()) ** 0.6
    inten = gaussian_filter(inten, sigma=max(0.6, 0.4 * scale))   # soft edges
    inten += 0.03 + 0.02 * rng.standard_normal(inten.shape)        # background + noise
    inten = np.clip(inten, 0, None)
    img = (inten / (inten.max() + 1e-9) * 60000).astype(np.uint16)
    return img, nlab


def main():
    ap = argparse.ArgumentParser(description="Render a synthetic Xenium folder for SegMap.")
    ap.add_argument("--gt-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--scale", type=float, default=3.0,
                    help="Pixels per synthetic length unit (default 3).")
    ap.add_argument("--format", choices=["parquet", "csv"], default="parquet",
                    help="transcripts file format (SegMap needs parquet).")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    tx = pd.read_csv(os.path.join(a.gt_dir, "transcripts.csv"))
    with open(os.path.join(a.gt_dir, "ground_truth_summary.json")) as f:
        summary = json.load(f)
    domain = summary["stats"]["domain"]
    W = int(np.ceil(domain[0] * a.scale)); H = int(np.ceil(domain[1] * a.scale))

    npz = np.load(os.path.join(a.gt_dir, "polygons.npz"))
    nuclei = {int(k.split("_")[1]): npz[k] for k in npz.files if k.startswith("nucleus_")}

    print(f"Rendering DAPI image {H}x{W} at scale={a.scale} ...")
    img, nlab = render_dapi(nuclei, H, W, a.scale, seed=a.seed)
    morph_path = os.path.join(a.out_dir, "morphology.ome.tif")
    tifffile.imwrite(morph_path, img)
    print(f"  wrote {morph_path}")

    # transcripts in Xenium schema, same pixel frame as the image
    out = pd.DataFrame({
        "transcript_id": np.arange(len(tx), dtype=np.int64),
        "x_location": tx["x"].to_numpy() * a.scale,
        "y_location": tx["y"].to_numpy() * a.scale,
        "z_location": 0.0,
        "feature_name": tx["gene"].astype(str).to_numpy(),
        "qv": 40.0,
        "overlaps_nucleus": 0,
    })
    tx_path = os.path.join(a.out_dir, f"transcripts.{a.format}")
    if a.format == "parquet":
        out.to_parquet(tx_path, index=False)
    else:
        out.to_csv(tx_path, index=False)
    print(f"  wrote {tx_path}  ({len(out):,} transcripts)")

    # recommended Cellpose diameter = median nucleus diameter in pixels
    diam = []
    for poly in nuclei.values():
        p = np.asarray(poly)
        diam.append(np.hypot(*(p.max(0) - p.min(0))) * a.scale * 0.5)  # ~mean extent
    med_diam = float(np.median(diam)) if diam else 30.0
    info = {"scale": a.scale, "image_H": H, "image_W": W, "domain": domain,
            "recommended_cellpose_diameter": round(med_diam, 1)}
    with open(os.path.join(a.out_dir, "synth_render_info.json"), "w") as f:
        json.dump(info, f, indent=2)

    print(f"\nRun SegMap with:  --cellpose_diameter {med_diam:.0f}  (median nucleus px)")
    print(f"  python -m segmap.run_segmap --xenium_dir {os.path.abspath(a.out_dir)} "
          f"--outdir {os.path.abspath(a.out_dir)}/segmap_out --no_tune_resolution "
          f"--cellpose_diameter {med_diam:.0f} ...")
    print(f"\nThen:  python segmap_to_pred.py --gt-dir {a.gt_dir} "
          f"--mask {os.path.abspath(a.out_dir)}/segmap_out/data/"
          f"full_mask_cellpose_fullres_tiled.tif --scale {a.scale} --out segmap_pred.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
