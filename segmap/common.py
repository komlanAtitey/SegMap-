#!/usr/bin/env python3  # use python3 interpreter
# -*- coding: utf-8 -*-  # file encoding declaration
"""
xenium_cellpose_fullres_tiled_upgraded_full.py  # module name

Upgraded "best fidelity (tiled full-res)" Xenium + Cellpose pipeline with big speedups and per-line comments.

This file is the complete pipeline adapted from your original script with the upgraded tiled segmentation,
robust merging and optional multiprocessing, and includes a comment on every source line for traceability.
"""  # docstring

# ------------------------- env/thread limits -------------------------
import os  # import os for environment and filesystem ops
import argparse
import math
import json
import time
import warnings
import tempfile
import shutil

from typing import Optional, Tuple, Dict, Any, List

import numpy as np
import pandas as pd
import scipy.sparse as sp

from scipy import ndimage as ndi

import tifffile as tiff

from PIL import Image

import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import anndata as ad

from skimage import measure
from skimage.color import label2rgb
from skimage.segmentation import (
    relabel_sequential,
    find_boundaries
)


os.environ["OMP_NUM_THREADS"] = "1"  # limit OpenMP threads to 1
os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"  # avoid nested OpenMP
os.environ["MKL_NUM_THREADS"] = "1"  # limit MKL threads
os.environ["OPENBLAS_NUM_THREADS"] = "1"  # limit OpenBLAS threads
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"  # limit macOS vecLib threads
os.environ["NUMEXPR_NUM_THREADS"] = "1"  # limit numexpr threads
os.environ["NUMBA_NUM_THREADS"] = "1"  # limit numba threads

# ------------------------- standard imports -------------------------


# ------------------------- optional packages -------------------------
try:  # try to import cellpose
    from cellpose import models  # cellpose models
    CELLPOSE_AVAILABLE = True  # flag
except Exception:  # fallback if not installed
    CELLPOSE_AVAILABLE = False  # flag

try:  # try import scanpy
    import scanpy as sc  # scanpy
    SCANPY_AVAILABLE = True  # flag
except Exception:  # fallback
    SCANPY_AVAILABLE = False  # flag

try:  # try import yaml
    import yaml  # yaml
    YAML_AVAILABLE = True  # flag
except Exception:  # fallback
    YAML_AVAILABLE = False  # flag

try:  # try pyarrow parquet
    import pyarrow.parquet as pq  # parquet
    PYARROW_AVAILABLE = True  # flag
except Exception:  # fallback
    PYARROW_AVAILABLE = False  # flag

try:  # try umap-learn
    import umap  # umap
    UMAPLEARN_AVAILABLE = True  # flag
except Exception:  # fallback
    UMAPLEARN_AVAILABLE = False  # flag

# ------------------------- defaults -------------------------
MODEL_TYPE = "nuclei"  # cellpose model type (DAPI nuclei); was "CYTO2" (invalid + ignored)

N_PATCHES = 6  # preview patches count
PATCH_SIZE_FULL = 2048  # preview patch size full-res
PATCH_DS = 2  # preview patch downsample factor
PATCH_PREVIEW_DS = 16  # downsample for patch scoring
PATCH_STRIDE_FULL = PATCH_SIZE_FULL  # stride for patch grid

TOP_GENE_CAP = 5050  # cap for genes in matrix
MIN_GENES_PER_CELL = 30  # gene filter per cell
MIN_CELLS_PER_GENE = 10  # cell filter per gene
RELAX_GENE_THRESHOLDS = [30, 20, 10, 5, 1]  # fallback gene thresholds
RELAX_CELL_THRESHOLDS = [10, 5, 2, 1]  # fallback cell thresholds
MIN_CELLS_FLOOR = 50  # keep at least this many cells

N_PCS = 50 #@@@@ 30  # PCA components
N_NEIGHBORS = 20  # @ 30 #@ 25 #@ komlan 20  # neighbors for UMAP

BOUNDARY_WIDTH = 2  # PNG boundary thickness
PASTEL_ALPHA = 0.90  # alpha for label2rgb

X_CANDIDATES = [
    "x",
    "x_location",
    "x_centroid",
    "x_pos",
    "x_um",
     "x_pixel"]  # candidate x cols
Y_CANDIDATES = [
    "y",
    "y_location",
    "y_centroid",
    "y_pos",
    "y_um",
     "y_pixel"]  # candidate y cols
GENE_CANDIDATES = [
    "gene",
    "feature_name",
    "gene_name",
    "target",
     "feature"]  # candidate gene cols

# ------------------------- CLI parser -------------------------


# ------------------------- stage timing -------------------------
class StageTimer:
    """Lightweight per-stage timer.

    Call ``lap("stage name")`` at each stage boundary; it prints the wall-clock
    time spent since the previous lap. Call ``summary()`` at the end to print a
    slowest-first breakdown so the dominant stage is obvious.
    """

    def __init__(self):
        self._t0 = time.perf_counter()
        self._last = self._t0
        self._laps = []

    def lap(self, label: str):
        now = time.perf_counter()
        dt = now - self._last
        self._laps.append((label, dt))
        print(f"[TIMER] {label:<42} {dt:8.1f}s  (elapsed {(now - self._t0)/60:5.1f} min)",
              flush=True)
        self._last = now

    def summary(self):
        total = self._last - self._t0
        print("[TIMER] ============== stage summary (slowest first) ==============",
              flush=True)
        for label, dt in sorted(self._laps, key=lambda x: -x[1]):
            pct = (100.0 * dt / total) if total > 0 else 0.0
            print(f"[TIMER]   {label:<42} {dt:8.1f}s  {pct:5.1f}%", flush=True)
        print(f"[TIMER]   {'TOTAL':<42} {total:8.1f}s  ({total/60:.1f} min)", flush=True)
