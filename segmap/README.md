# segmap_pkg

Modular refactor of the Xenium tiled-Cellpose segmentation + clustering pipeline.
Behavior is unchanged from the single-file `segmap_atitey_upgraded.py`; the code is
split into logical modules with no circular imports.

## Run

```bash
python run_segmap.py [same CLI args as before]   # e.g. --input ... --outdir ...
```

`run_segmap.py` simply calls `segmap_pkg.pipeline.main()`. All previous flags
(`--resolution`, `--tune_resolution/--no_tune_resolution`, `--spatial_weight`,
`--cellpose_diameter`, `--tile_batch`, …) are unchanged.

## Layout

| Module             | Responsibility |
|--------------------|----------------|
| `common.py`        | All third-party imports, thread-limit env vars, optional-dependency guards (`*_AVAILABLE` flags) and tunable constants. Every other module does `from .common import *`. |
| `utils.py`         | Generic helpers: `pick_first_present`, `normalize01`, `downsample_mean`. |
| `io_ome.py`        | OME-TIFF plane loading, nucleus-plane scoring, patch/mosaic previews. |
| `segmentation.py`  | `CellposeRunner` (+ batched eval), multiprocessing workers, tiled segmentation, mask PNG export. |
| `assignment.py`    | Transcript→cell assignment and cell-by-gene matrix/CSV export. |
| `preprocessing.py` | AnnData QC filtering and PCA/HVG embedding. |
| `clustering.py`    | Subsampled silhouette helpers and Leiden (igraph) with resolution tuning. |
| `markers.py`       | Fast fold-change markers, Scanpy Wilcoxon markers, auto-annotation. |
| `refinement.py`    | SegMap Bayesian spatial refinement and stage snapshots. |
| `visualization.py` | UMAP / spatial plots and cluster CSV exports. |
| `cli.py`           | `parse_args`. |
| `pipeline.py`      | `main()` orchestration. |

## Internal dependencies (acyclic)

```
common  <-  (everything)
utils   <-  io_ome, segmentation
clustering, markers, assignment  <-  refinement
all modules  <-  pipeline
```

## Note

`pipeline.py` (~900 lines) still holds the whole `main()`. A natural next step is to
break it into stage functions (load → segment → assign → cluster → refine → plot),
but that is a behavioral refactor and was left out of this purely structural split.

## Output layout

`<xenium_dir>/cellpose_pipeline_outputs/` is split into two groups:

```
cellpose_pipeline_outputs/
├── figures/                      # group 1: all figures (PNG)
│   ├── umap_clusters.png
│   ├── umap_celltypes_unsupervised.png
│   ├── spatial_clusters.png, spatial_celltypes_unsupervised.png
│   ├── cellpose_polygons_mosaic_*.png  + patch previews
│   └── segmentation overlays
└── data/                         # group 2: all numerical data
    ├── xenium_cellpose_fullres_tiled_final.h5ad   (canonical AnnData)
    ├── full_mask_cellpose_fullres_tiled.tif       (label mask)
    ├── cell_by_gene_for_umap*.csv
    ├── cluster_markers_all_*.csv, cluster_labels_unsupervised.json
    ├── spatial_clusters*.csv, umap_clusters*.csv, posterior_entropy_final.csv
    ├── transcript_assignment_summary.json, pipeline_summary.json
    └── seurat/                                     (Seurat object bundle)
        ├── matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz
        ├── metadata.csv, pca.csv, umap.csv, umap_initial.csv, spatial.csv
        └── make_seurat.R
```

### Building the Seurat object

A native Seurat `.rds` is an R object and can't be written from Python, so the
pipeline writes a Seurat-loadable bundle. Build the object once in R:

```bash
cd <xenium_dir>/cellpose_pipeline_outputs/data/seurat
Rscript make_seurat.R        # needs install.packages(c("Seurat","Matrix"))
# -> seurat_object.rds  (counts + all metadata + pca/umap/spatial reductions)
```

## Cell types vs. spatial domains

The pipeline now produces two distinct labels instead of letting the spatial
refinement overwrite cell identity:

- **`cluster` / `cluster_final`** — transcriptomic cell types. These are the
  `cluster_initial` Leiden clusters, then **spatially denoised**: a low-confidence
  cell is reassigned to its spatial-neighborhood majority only when a strong
  majority of its neighbors agree. Confident/interior cells are never changed, so
  identity is preserved and only isolated misassignments are cleaned up. The final
  UMAP (`X_umap_final`) is the transcriptomic embedding, so it is as clean as the
  initial one.
- **`spatial_domain`** — niche / spatial-domain label from the Leiden clustering on
  the hybrid `X_bayes` space. Kept as a separate column (and `figures/spatial_domains.png`),
  not used as the cell type.

Relevant flags:

- `--spatial_denoise_iters` (default 2; `0` disables denoising)
- `--spatial_denoise_agree` (default 0.65) — neighbor-agreement threshold to reassign a low-confidence cell
- `--spatial_weight` (default 0.0) — weight of explicit spatial/density terms in the `spatial_domain` clustering space

## Benchmark-oriented improvements (opt-in; defaults preserve prior behavior)

These target marker strength and spatial coherence. A/B each against the
benchmark rather than enabling all at once:

- `--normalization pearson` — analytic Pearson residuals instead of CP10k+log1p+scale
  (better separability / marker contrast for sparse count panels). Also fixes the
  prior HVG bug (seurat_v3 flavor was run on log-normalized data).
- `--auto_n_pcs` — choose the PC count from the variance-ratio elbow (drops noise PCs).
- `--subtract_background` — subtract ambient signal estimated from Xenium
  negative-control probes/codewords before clustering, then drop control features.
- `--remove_doublets` — drop predicted doublets (Scrublet) before clustering.

Not yet implemented (largest lever, needs your morphology channels): whole-cell
(not nucleus-only) segmentation using non-DAPI morphology stains.
