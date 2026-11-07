"""Transcript-to-cell assignment and cell-by-gene matrix export."""
from .common import *  # noqa: F401,F403


def save_cell_by_gene_csv(
    adata,
    out_csv_path,
    use_layer: Optional[str] = None,
    index: bool = True,
     compress: bool = False):
    """
    Save cell-by-gene expression matrix suitable for UMAP input.

    Parameters
    ----------
    adata : anndata.AnnData or pandas.DataFrame
        If AnnData, will use adata.X (or adata.layers[use_layer] if provided).
        If DataFrame, will be used directly.
    out_csv_path : str
        Destination CSV path (or .csv.gz if compress=True).
    use_layer : str or None
        If provided and adata has .layers, use that layer.
    index : bool
        Whether to write the row index (cell IDs).
    compress : bool
        If True, write gzipped CSV (adds .gz if not present).
    """
    import pandas as pd
    # @@@@@@@@@@ komlan import numpy as np
    from scipy import sparse
    from pathlib import Path

    out_p = Path(out_csv_path)
    if compress and not out_p.suffix == ".gz":
        out_p = out_p.with_suffix(out_p.suffix + ".gz")

    # Build DataFrame from AnnData or accept DataFrame directly
    if hasattr(
    adata,
    "to_df") or (
        hasattr(
            adata,
            "X") and hasattr(
                adata,
                 "obs")):
        # choose layer if requested
        if use_layer is not None and hasattr(
    adata, "layers") and use_layer in adata.layers:
            mat = adata.layers[use_layer]
        else:
            try:
                # prefer a labeled DataFrame if available
                mat = adata.to_df()
            except Exception:
                mat = adata.X

        if isinstance(mat, pd.DataFrame):
            df = mat.copy()
        else:
            # handle sparse/dense numpy arrays
            if sparse.issparse(mat):
                mat = mat.toarray()
            # try to get row/col names
            try:
                obs_idx = getattr(
    adata, "obs_names", None) or getattr(
        adata, "obs", None).index
            except Exception:
                obs_idx = None
            try:
                var_idx = getattr(
    adata, "var_names", None) or getattr(
        adata, "var", None).index
            except Exception:
                var_idx = None
            df = pd.DataFrame(mat, index=obs_idx, columns=var_idx)
    elif isinstance(adata, pd.DataFrame):
        df = adata.copy()
    else:
        raise ValueError("adata must be an AnnData or a pandas DataFrame")

    # Keep only numeric columns (UMAP expects numeric features)
    numeric_df = df.select_dtypes(include=[np.number]).fillna(0)

    # If there are no column names (rare), create generic gene names
    if numeric_df.columns.isnull().any():
        numeric_df.columns = [f"gene_{i}" for i in range(numeric_df.shape[1])]

    # Ensure directory exists
    out_p.parent.mkdir(parents=True, exist_ok=True)

    # Write CSV; compress if requested
    if compress:
        numeric_df.to_csv(out_p, index=index, compression="gzip")
    else:
        numeric_df.to_csv(out_p, index=index)

    print(f"Saved cell-by-gene CSV for UMAP: {str(out_p)}", flush=True)


def parquet_schema_columns(path: str) -> Optional[List[str]]:
    if not PYARROW_AVAILABLE:  # no pyarrow
        return None  # return none
    pf = pq.ParquetFile(path)  # open parquet
    return list(pf.schema.names)  # return column names


def assign_transcripts_to_cells_fullres(mask_full: np.ndarray, tx: pd.DataFrame, x_col: str, y_col: str) -> Tuple[pd.DataFrame, Dict[str, Any]]:  # assign
    Hm, Wm = mask_full.shape  # dims
    x = tx[x_col].to_numpy(dtype=np.float64, copy=False)  # x array
    y = tx[y_col].to_numpy(dtype=np.float64, copy=False)  # y array

    def _try_map(xp: np.ndarray, yp: np.ndarray, label: str) -> Dict[str, Any]:  # mapping attempt
        xi = np.rint(xp).astype(np.int32, copy=False)  # round x
        yi = np.rint(yp).astype(np.int32, copy=False)  # round y
        valid = (xi >= 0) & (yi >= 0) & (xi < Wm) & (yi < Hm)  # bounds
        cell_id = np.zeros(xi.shape[0], dtype=np.int32)  # init
        if np.any(valid):  # any valid
            cell_id[valid] = mask_full[yi[valid], xi[valid]].astype(np.int32, copy=False)  # lookup
        n_valid = int(valid.sum())  # count
        n_assigned = int((cell_id > 0).sum())  # assigned
        n_cells = int(np.unique(cell_id[cell_id > 0]).size) if n_assigned > 0 else 0  # unique cells
        return {"label": label, "cell_id": cell_id, "n_valid": n_valid, "n_assigned": n_assigned, "n_cells": n_cells}  # result

    x_min, x_max = float(np.nanmin(x)), float(np.nanmax(x))  # ranges
    y_min, y_max = float(np.nanmin(y)), float(np.nanmax(y))  # ranges
    print(f"[TX] Raw coord ranges: x=[{x_min:.3f},{x_max:.3f}] y=[{y_min:.3f},{y_max:.3f}]", flush=True)  # log
    print(f"[MASK] full-res shape H={Hm} W={Wm}", flush=True)  # log

    candidates: List[Dict[str, Any]] = []  # candidates
    candidates.append(_try_map(x, y, "raw_pixels"))  # raw
    candidates.append(_try_map(x - x_min, y - y_min, "shift_to_zero"))  # shift
    candidates.append(_try_map(y, x, "swap_xy_raw"))  # swap
    candidates.append(_try_map(y - y_min, x - x_min, "swap_xy_shift_to_zero"))  # swap+shift

    xr = float(x_max - x_min)  # x range
    yr = float(y_max - y_min)  # y range
    if (xr > 0) and (yr > 0):  # guard
        sx = float(Wm / xr)  # scale guess x
        sy = float(Hm / yr)  # scale guess y
        if (sx > 1.2) and (sy > 1.2):  # likely microns -> pixels
            candidates.append(_try_map((x - x_min) * sx, (y - y_min) * sy, f"shift+scale(sx={sx:.3f},sy={sy:.3f})"))  # add scaled
            candidates.append(_try_map((y - y_min) * sy, (x - x_min) * sx, f"swap+shift+scale(sx={sx:.3f},sy={sy:.3f})"))  # add scaled swap

    best = max(candidates, key=lambda d: (d["n_assigned"], d["n_valid"], d["n_cells"]))  # choose best
    print(f"[ASSIGN] Best mapping='{best['label']}' n_valid={best['n_valid']} n_assigned={best['n_assigned']} unique_cells={best['n_cells']}", flush=True)  # log

    tx2 = tx.copy()  # copy
    tx2["cell_id"] = best["cell_id"]  # assign
    tx_assigned = tx2[tx2["cell_id"] > 0].copy()  # keep assigned

    summary = {k: best[k] for k in ("label", "n_valid", "n_assigned", "n_cells")}  # summary
    summary.update({"mask_H": int(Hm), "mask_W": int(Wm), "raw_x_min": x_min, "raw_x_max": x_max, "raw_y_min": y_min, "raw_y_max": y_max})  # add
    return tx_assigned, summary  # return
