"""Export the final AnnData as a Seurat-loadable bundle + an R script.

A native Seurat ``.rds`` is an R (S4) object and cannot be written reliably from
Python, so this writes the numerical data in a form Seurat reads directly:

    seurat/
      matrix.mtx.gz        counts, features x cells (10x convention)
      barcodes.tsv.gz      cell ids
      features.tsv.gz      gene ids (id, name, type)
      metadata.csv         per-cell obs (clusters, cell types, QC, coords, ...)
      pca.csv              PCA embedding (if present)
      umap.csv             final UMAP embedding (if present)
      umap_initial.csv     transcriptomic UMAP (if present)
      spatial.csv          centroid coordinates (if present)
      make_seurat.R        run `Rscript make_seurat.R` to build seurat_object.rds

Running the R script once yields a real Seurat object containing the counts,
all per-cell metadata, and the embeddings as DimReduc objects.
"""
from .common import *  # noqa: F401,F403

import gzip
import scipy.io


def _write_embedding(obsm_key: str, adata, out_dir: str, fname: str, prefix: str) -> bool:
    if obsm_key not in adata.obsm:
        return False
    Y = np.asarray(adata.obsm[obsm_key])
    cols = [f"{prefix}{i+1}" for i in range(Y.shape[1])]
    df = pd.DataFrame(Y, columns=cols)
    df.insert(0, "cell_id", np.asarray(adata.obs_names).astype(str))
    df.to_csv(os.path.join(out_dir, fname), index=False)
    return True


def export_seurat_bundle(adata, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    # ---- counts: Seurat/10x store features x cells ----
    X = adata.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    counts = X.T.tocoo()  # genes (rows) x cells (cols)
    counts = counts.astype(np.int64) if np.allclose(counts.data, np.round(counts.data)) else counts
    with gzip.open(os.path.join(out_dir, "matrix.mtx.gz"), "wb") as fh:
        scipy.io.mmwrite(fh, counts, field="integer"
                         if counts.dtype.kind in "iu" else "real")

    # ---- barcodes (cell ids) ----
    pd.Series(np.asarray(adata.obs_names).astype(str)).to_csv(
        os.path.join(out_dir, "barcodes.tsv.gz"),
        index=False, header=False, compression="gzip")

    # ---- features (gene id, name, type) — 10x 3-column format ----
    genes = np.asarray(adata.var_names).astype(str)
    pd.DataFrame({"id": genes, "name": genes,
                  "type": "Gene Expression"}).to_csv(
        os.path.join(out_dir, "features.tsv.gz"),
        sep="\t", index=False, header=False, compression="gzip")

    # ---- per-cell metadata ----
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", np.asarray(adata.obs_names).astype(str))
    obs.to_csv(os.path.join(out_dir, "metadata.csv"), index=False)

    # ---- embeddings ----
    _write_embedding("X_pca_fast", adata, out_dir, "pca.csv", "PC_")
    _write_embedding("X_umap_final", adata, out_dir, "umap.csv", "UMAP_")
    _write_embedding("X_umap_fast", adata, out_dir, "umap_initial.csv", "UMAPi_")
    # spatial coordinates from centroid columns, if present
    sp_cols = [c for c in ("x_centroid_px", "y_centroid_px") if c in adata.obs.columns]
    if len(sp_cols) == 2:
        sdf = adata.obs[sp_cols].copy()
        sdf.columns = ["Spatial_1", "Spatial_2"]
        sdf.insert(0, "cell_id", np.asarray(adata.obs_names).astype(str))
        sdf.to_csv(os.path.join(out_dir, "spatial.csv"), index=False)

    _write_r_script(out_dir)
    print(f"[SEURAT] Wrote Seurat bundle to {out_dir} "
          f"(build with: Rscript {os.path.join(out_dir, 'make_seurat.R')})",
          flush=True)


_R_SCRIPT = r'''#!/usr/bin/env Rscript
# Build a Seurat object from the exported bundle.
#   Rscript make_seurat.R
# Requires: install.packages(c("Seurat","Matrix"))
suppressMessages({library(Seurat); library(Matrix)})

# work in this script's own directory
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("--file=", "", grep("--file=", args, value = TRUE)))
if (length(here) == 1 && nzchar(here)) setwd(here)

counts <- ReadMtx(mtx = "matrix.mtx.gz",
                  cells = "barcodes.tsv.gz",
                  features = "features.tsv.gz",
                  feature.column = 1)

meta <- read.csv("metadata.csv", row.names = 1, check.names = FALSE)
meta <- meta[colnames(counts), , drop = FALSE]

obj <- CreateSeuratObject(counts = counts, meta.data = meta)

add_reduction <- function(obj, file, key) {
  if (!file.exists(file)) return(obj)
  emb <- read.csv(file, row.names = 1, check.names = FALSE)
  emb <- as.matrix(emb[colnames(obj), , drop = FALSE])
  obj[[key]] <- CreateDimReducObject(
    embeddings = emb, key = paste0(toupper(key), "_"),
    assay = DefaultAssay(obj))
  obj
}
obj <- add_reduction(obj, "pca.csv", "pca")
obj <- add_reduction(obj, "umap.csv", "umap")
obj <- add_reduction(obj, "umap_initial.csv", "umapinitial")
obj <- add_reduction(obj, "spatial.csv", "spatial")

saveRDS(obj, "seurat_object.rds")
cat("Saved seurat_object.rds:", ncol(obj), "cells x", nrow(obj), "features\n")
'''


def _write_r_script(out_dir: str) -> None:
    with open(os.path.join(out_dir, "make_seurat.R"), "w") as f:
        f.write(_R_SCRIPT)
