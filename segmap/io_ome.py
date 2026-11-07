"""OME-TIFF plane loading, nucleus-plane scoring, and patch/mosaic previews."""
from .common import *  # noqa: F401,F403
from .utils import normalize01


def score_plane_for_nuclei(plane: np.ndarray, ds: int = 16) -> float:
    p = plane[::ds, ::ds].astype(np.float32, copy=False)  # downsample
    p = ndi.gaussian_filter(p, 1.0)  # gaussian smooth
    gx = ndi.sobel(p, axis=1)  # sobel x
    gy = ndi.sobel(p, axis=0)  # sobel y
    g = np.hypot(gx, gy)  # gradient magnitude
    return float(g.var())  # return variance


def _axes_index_map(axes: str) -> Dict[str, int]:  # map axis char to index
    return {a: i for i, a in enumerate(list(axes))}  # dict


def load_one_plane_ome(ome_tif_path, auto_pick=True, score_ds=16, zproject="none"):
    """
import shutil  # shutil for cleanup
import tempfile  # import tempfile for memmap
from sklearn.cluster import KMeans  # KMeans clustering
from sklearn.decomposition import PCA  # PCA
import matplotlib.pyplot as plt  # plotting
import anndata as ad  # AnnData container
import scipy.sparse as sp  # sparse matrices
from PIL import Image  # image saving
from skimage.segmentation import relabel_sequential  # relabel function
from skimage.segmentation import find_boundaries  # find mask boundaries
from skimage.color import label2rgb  # colored label plotting
from skimage import measure  # regionprops and measure tools
from scipy import ndimage as ndi  # ndimage operations
import tifffile as tiff  # tiff image IO
import pandas as pd  # dataframes
import numpy as np  # numerical arrays
from typing import Optional, Tuple, Dict, Any, List  # typing aids
import warnings  # warnings control
import json  # json I/O
import math  # math helpers
import time  # timing utilities
import argparse  # CLI argument parsing
    Load one best DAPI plane from an OME-TIFF.  # doc
    Supports axes like:  # doc
      - ZYX  # doc
      - CZYX  # doc
      - TCZYX  # doc
      - YX  # doc
    Returns:  # doc
      dapi_2d (float32), chosen_plane (int), axes_string (str)  # doc
    """  # docstring

    import tifffile  # local import
    # @@@@@@@@@@ komlan import numpy as np  # local import

    with tifffile.TiffFile(ome_tif_path) as tif:  # open OME-TIFF
        series = tif.series[0]  # take first series
        axes = series.axes  # axes string
        shape = series.shape  # array shape

        print(f"OME axes: {axes} shape: {shape}")  # log axes/shape

        arr = series.asarray()  # load full array

    arr = np.asarray(arr)  # ensure numpy array

    ax = {a: i for i, a in enumerate(axes)}  # axis index map

    if "Y" not in ax or "X" not in ax:  # require Y/X
        raise ValueError(
    f"OME-TIFF missing Y/X axes. Got axes={axes}")  # error

    if "C" in ax:  # if channel axis exists
        cdim = arr.shape[ax["C"]]  # channel dim
        c_idx = 0  # Xenium usually: DAPI is channel 0
    else:  # no channel axis
        c_idx = None  # no channel selection

    if "T" in ax:  # if time axis exists
        t_idx = 0  # first timepoint
    else:  # no time axis
        t_idx = None  # no time selection

    if "Z" in ax:  # if Z axis exists
        zdim = arr.shape[ax["Z"]]  # number of planes
    else:  # no Z axis
        zdim = 1  # single plane

    def get_plane(z):  # helper: slice out a 2D plane
        sl = [slice(None)] * arr.ndim  # start with all slices

        if t_idx is not None:  # if time axis
            sl[ax["T"]] = t_idx  # set time index
        if c_idx is not None:  # if channel axis
            sl[ax["C"]] = c_idx  # set channel index
        if "Z" in ax:  # if Z axis
            sl[ax["Z"]] = z  # set z plane index

        plane = arr[tuple(sl)]  # slice array

        plane = np.squeeze(plane)  # drop singleton dims

        if plane.ndim != 2:  # ensure 2D
            raise ValueError(
    f"Expected 2D plane, got shape={
        plane.shape}, axes={axes}")  # error

        return plane  # return plane

    if zdim == 1:  # single plane
        dapi = get_plane(0).astype(np.float32)  # cast to float32
        return dapi, 0, axes  # return

    # ---- Optional Z-projection across DAPI focal planes (this image is a DAPI
    # z-stack). Collapsing focal planes gives Cellpose a sharper, more complete
    # nucleus image than any single plane, recovering nuclei that are only in
    # focus at some depths. Streamed plane-by-plane to bound memory.
    zp = str(zproject).lower()
    if zp in ("max", "mean") and zdim > 1:
        acc = get_plane(0).astype(np.float32)
        if zp == "max":
            for z in range(1, zdim):
                np.maximum(acc, get_plane(z).astype(np.float32, copy=False), out=acc)
        else:  # mean
            for z in range(1, zdim):
                acc += get_plane(z)
            acc /= float(zdim)
        print(f"[DAPI] z-projection '{zp}' over {zdim} planes -> shape {acc.shape}",
              flush=True)
        return acc.astype(np.float32, copy=False), -1, axes  # -1 = projected (no single plane)

    if auto_pick:  # auto select best focus
        best_z = 0  # initialize best z
        best_score = -np.inf  # initialize best score

        for z in range(zdim):  # loop over planes
            p = get_plane(z)  # get plane

            if score_ds > 1:  # downsample for speed
                p_small = p[::score_ds, ::score_ds]  # subsample
            else:  # no downsample
                p_small = p  # use original

            gy, gx = np.gradient(p_small.astype(np.float32))  # gradient
            score = float(np.mean(gx * gx + gy * gy))  # focus score

            if score > best_score:  # if better score
                best_score = score  # update score
                best_z = z  # update best z

        chosen_plane = best_z  # chosen plane
    else:  # not auto pick
        chosen_plane = zdim // 2  # middle plane

    dapi = get_plane(chosen_plane).astype(np.float32)  # extract chosen plane
    return dapi, chosen_plane, axes  # return


def pick_top_patches_by_signal(dapi_full: np.ndarray,
    patch_size: int,
    n_patches: int,
    stride: int,
    preview_ds: int = 16) -> List[Tuple[int,
     int]]:
    d = dapi_full[::preview_ds, ::preview_ds]  # downsample
    d = normalize01(d)  # normalize
    hp = d - ndi.gaussian_filter(d, 6.0)  # high-pass
    hp = np.clip(hp, 0.0, 1.0)  # clip
    H, W = dapi_full.shape  # full dims
    coords: List[Tuple[int, int]] = []  # coord list
    scores: List[float] = []  # scores
    for y0 in range(0, max(1, H - patch_size), stride):  # grid y
        for x0 in range(0, max(1, W - patch_size), stride):  # grid x
            ys, xs = y0 // preview_ds, x0 // preview_ds  # preview start
            ye = min(
    (y0 + patch_size) // preview_ds,
     hp.shape[0])  # preview end y
            xe = min(
    (x0 + patch_size) // preview_ds,
     hp.shape[1])  # preview end x
            win = hp[ys:ye, xs:xe]  # window
            if win.size == 0:  # skip empty
                continue  # continue
            coords.append((x0, y0))  # append coord
            scores.append(float(win.mean()))  # append score
    if len(scores) == 0:  # no scores
        return [(0, 0)] * int(n_patches)  # default
    idx = np.argsort(np.asarray(scores))[::-1][: int(n_patches)]  # top indices
    return [coords[i] for i in idx]  # return top coords


def save_patch_polygon_image(label_mask: np.ndarray, out_png: str, boundary_width: int = 2, pastel_alpha: float = 0.9) -> None:  # save overlay
    label_mask = label_mask.astype(np.int32, copy=False)  # ensure int
    if int(label_mask.max()) == 0:  # empty mask
        blank = np.zeros(
    (label_mask.shape[0],
    label_mask.shape[1],
    3),
     dtype=np.uint8)  # blank image
        Image.fromarray(blank).save(out_png)  # save
        print("WARNING: empty mask saved:", out_png, flush=True)  # warn
        return  # done
    rgb = (
    label2rgb(
        label_mask,
        bg_label=0,
        alpha=float(pastel_alpha)) *
        255).astype(
            np.uint8)  # color labels
    b = find_boundaries(label_mask, mode="outer")  # find boundaries
    if int(boundary_width) > 1:  # thicken
        b = ndi.binary_dilation(
    b, iterations=int(boundary_width) - 1)  # dilate
    rgb[b] = (0, 0, 0)  # black boundary
    Image.fromarray(rgb).save(out_png, format="PNG")  # save png
    print("Saved patch:", out_png, "cells:", int(
        label_mask.max()), flush=True)  # log


def save_mosaic(png_paths: List[str], out_png: str, cols: int = 3, padding: int = 12, bg: Tuple[int, int, int] = (255, 255, 255)) -> None:  # mosaic
    imgs = [Image.open(p) for p in png_paths]  # load images
    if len(imgs) == 0:  # none
        raise RuntimeError("No patch PNGs to mosaic.")  # error
    tile_w = max(im.width for im in imgs)  # tile width
    tile_h = max(im.height for im in imgs)  # tile height
    cols = int(min(int(cols), len(imgs)))  # cols
    rows = int(math.ceil(len(imgs) / cols))  # rows
    canvas_w = cols * tile_w + (cols + 1) * int(padding)  # canvas width
    canvas_h = rows * tile_h + (rows + 1) * int(padding)  # canvas height
    canvas = Image.new("RGB", (canvas_w, canvas_h), color=bg)  # new canvas
    for i, im in enumerate(imgs):  # iterate
        r = i // cols  # row
        c = i % cols  # col
        x = int(padding) + c * (tile_w + int(padding))  # x pos
        y = int(padding) + r * (tile_h + int(padding))  # y pos
        ox = x + (tile_w - im.width) // 2  # offset x
        oy = y + (tile_h - im.height) // 2  # offset y
        canvas.paste(im, (ox, oy))  # paste
    canvas.save(out_png, format="PNG")  # save mosaic
    print("Saved mosaic:", out_png, flush=True)  # log
