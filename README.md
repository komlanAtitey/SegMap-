# SegMap

**SegMap** is an integrated computational framework for **nucleus-guided cell reconstruction, transcript assignment, spatial marker identification, and cell–cell network analysis** in imaging-based spatial transcriptomics. Starting from Xenium morphology images and transcript coordinates, SegMap reconstructs cells using physically observed nuclei as spatial anchors, expands nuclear boundaries to establish cellular territories, resolves overlaps, and assigns transcripts to reconstructed cells. It subsequently constructs spatial neighborhood graphs and identifies spatially informative marker genes through Bayesian spatial modeling. Importantly, SegMap incorporates uncertainty-aware iterative refinement, using transcript-assignment uncertainty and spatial marker evidence to improve cell boundaries and cell-state annotations.

Unlike transcript-density-based segmentation methods, SegMap generates **compact, non-overlapping cell territories** while preserving acellular regions, such as vascular lumina and alveolar air spaces. By integrating cellular reconstruction, transcript assignment, spatially informed clustering, Bayesian marker discovery, and uncertainty-guided refinement within a unified framework, SegMap provides a robust foundation for characterizing cellular heterogeneity, spatial organization, and cell–cell relationships in complex tissues.

---

# Overview

## SegMap Workflow

```text
Xenium Data
    │
    ▼
Nucleus Segmentation (Cellpose)
    │
    ▼
Cell Reconstruction
    │
    ▼
Transcript Assignment
    │
    ▼
Cell–Cell Spatial Network
    │
    ▼
Bayesian Spatial Marker Analysis
    │
    ▼
Spatial Clustering & Visualization
```

The complete workflow consists of five stages:

### 1. Nucleus Segmentation

Nuclei are segmented directly from Xenium DAPI morphology images using the Cellpose nuclei model.

**Input**

```text
morphology_focus/
morphology.ome.tif
```

**Output**

```text
nucleus_mask.tif
nucleus_polygons.csv
```

---

### 2. Cell Reconstruction

Each segmented nucleus is expanded by a user-defined radius to generate a reconstructed cell territory. Cell territories are mutually exclusive and preserve background regions.

**Output**

```text
reconstructed_cells.tif
cell_boundaries.csv
cell_morphology.csv
```

---

### 3. Transcript Assignment

Decoded transcripts are assigned to reconstructed cells according to their spatial coordinates.

**Output**

```text
cell_gene_matrix.h5ad
cell_gene_matrix.csv
```

---

### 4. Cell–Cell Spatial Network Construction

Cell centroids are used to construct a k-nearest-neighbor graph representing local tissue architecture.

**Output**

```text
spatial_graph.edgelist
spatial_neighbors.csv
```
---

### 5. Uncertainty-Aware Iterative Refinement

SegMap incorporates an iterative refinement procedure that is activated when transcript-assignment uncertainty exceeds a predefined threshold. The procedure uses assignment uncertainty and spatial information to identify potentially inaccurate cell boundaries and ambiguous cell-state assignments. Cell territories and transcript assignments are subsequently refined while preserving non-overlapping spatial organization and biologically plausible cellular geometry. The iterative procedure terminates when the uncertainty criterion is satisfied or the predefined iteration limit is reached.

**Output**

```text
refined_cells.tif
refined_cell_boundaries.csv
refined_cell_morphology.csv
refined_cell_gene_matrix.h5ad
refined_cell_gene_matrix.csv
refinement_summary.json
```

---

### 6. Bayesian Spatial Marker Analysis

SegMap identifies genes exhibiting statistically significant spatial organization and cluster-specific enrichment.

For each gene, SegMap estimates:

- Posterior probability of spatial relevance
- Fraction of Spatial Variance (FSV)
- Bayes Factor
- Spatial effect size
- Cluster-specific spatial markers

**Output**

```text
spatial_markers.csv
cluster_markers.csv
marker_posteriors.csv
```
---

# Repository Structure

```text
SegMap/
│
├── README.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── setup.py
├── pyproject.toml
│
├── segmap/
│   ├── run_segmap.py
│   ├── segmentation/
│   ├── reconstruction/
│   ├── transcript_assignment/
│   ├── spatial_graph/
│   ├── bayesian_marker_model/
│   ├── clustering/
│   ├── visualization/
│   └── utils/
│
├── data/
│   ├── raw/
│   ├── processed/
│   └── benchmark/
│
├── simulations/
│   ├── synthetic_generation/
│   ├── benchmark_framework/
│   └── evaluation_metrics/
│
├── benchmarking/
│   ├── segmap/
│   ├── scanpy/
│   ├── seurat/
│   ├── baysor/
│   └── comparison_results/
│
├── analysis/
│   ├── R/
│   └── Python/
│
├── figures/
│   ├── manuscript/
│   ├── supplementary/
│   └── source_data/
│
├── results/
│   ├── datasets/
│   ├── simulations/
│   └── benchmarks/
│
├── notebooks/
│
└── docs/
    ├── methods/
    ├── tutorials/
    └── examples/
```

---

# Installation

## Clone the Repository

```bash
git clone https://github.com/komlanAtitey/SegMap-.git
cd SegMap-
```

## Create a Virtual Environment

### Linux / macOS

```bash
python -m venv segmap_env
source segmap_env/bin/activate
```

### Windows

```bash
python -m venv segmap_env
segmap_env\Scripts\activate
```

## Install Dependencies

```bash
pip install -r requirements.txt
```

---

# Required Input

SegMap expects a standard Xenium output directory:

```text
XeniumData/
└── Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs/
    ├── cells.csv.gz
    ├── transcripts.csv.gz
    ├── morphology_focus/
    ├── morphology.ome.tif
    ├── experiment.xenium
    └── ...
```

---

# Running SegMap

```bash
python -m segmap.run_segmap \
  --xenium_dir XeniumData/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs \
  --reuse_mask auto \
  --gpu \
  --no_tune_resolution \
  --resolution 0.3 \
  --n_pcs 20 \
  --expand_pixels 80 \
  --min_genes_per_cell 10 \
  --min_counts_per_cell 25 \
  --min_cells_per_gene 3 \
  --top_gene_cap 5050 \
  --cellpose_diameter 36 \
  --umap_init paga \
  --umap_min_dist 0.1
```

---

# SegMap Pipeline Outputs

Running SegMap generates two main output directories: `figures/` and `data/`.

```text
cellpose_pipeline_outputs/
│
├── figures/
│   ├── cellpose_polygons_mosaic_6patches.png
│   ├── cellpose_polygons_patch01.png
│   ├── cellpose_polygons_patch02.png
│   ├── cellpose_polygons_patch03.png
│   ├── cellpose_polygons_patch04.png
│   ├── cellpose_polygons_patch05.png
│   ├── cellpose_polygons_patch06.png
│   ├── segmentation_boundaries_fullres.png
│   ├── segmentation_colored_fullres.png
│   ├── spatial_celltypes_unsupervised.png
│   ├── spatial_clusters.png
│   ├── spatial_domains.png
│   ├── umap_celltypes_unsupervised.png
│   └── umap_clusters.png
│
└── data/
    ├── cellpose_cells_table_fullres.csv
    ├── spatial_clusters_final.csv
    ├── umap_clusters_final.csv
    ├── cluster_markers_all_fast_final.csv
    ├── cluster_markers_all_cells_final_top20.csv
    ├── posterior_entropy_final.csv
    ├── Benchmark_bayes_markers.csv
    ├── cell_by_gene_for_umap_final.csv
    ├── cluster_separability.csv
    ├── cluster_labels_unsupervised.json
    ├── pipeline_summary.json
    ├── overlap_resolution_qc.json
    └── ...
```

#### Figures

The `figures/` directory contains visual outputs generated by the SegMap pipeline, including:

- **Cell segmentation:** Cellpose polygon visualizations across six image patches, segmentation boundaries, and full-resolution colored segmentation.
- **Spatial visualization:** Spatial distribution of unsupervised cell-type assignments, clusters, and spatial domains.
- **Dimensionality reduction:** UMAP visualizations of unsupervised cell-type assignments and clusters.

#### Data

The `data/` directory contains quantitative results, intermediate datasets, and pipeline reports, including:

- **Cell segmentation:** Full-resolution Cellpose cell table.
- **Clustering:** Spatial coordinates, UMAP coordinates, cluster assignments, and unsupervised cluster labels.
- **Marker identification:** Cluster-specific marker genes and Bayesian marker analysis results.
- **Uncertainty quantification:** Posterior entropy estimates.
- **Cluster evaluation:** Cluster separability metrics.
- **Quality control and reproducibility:** Pipeline execution summary and overlap-resolution quality-control metrics.

Additional data files generated during pipeline execution may also be stored in this directory.

---
# Key Parameters

| Parameter | Description |
|------------|------------|
| `--xenium_dir` | Path to Xenium output directory |
| `--gpu` | Enable GPU acceleration |
| `--reuse_mask auto` | Reuse existing segmentation if available |
| `--resolution` | Working image resolution (µm/pixel) |
| `--cellpose_diameter` | Expected nuclear diameter for Cellpose |
| `--expand_pixels` | Cell reconstruction radius around nuclei |
| `--n_pcs` | Number of principal components |
| `--min_genes_per_cell` | Minimum genes required per cell |
| `--min_counts_per_cell` | Minimum transcript counts per cell |
| `--min_cells_per_gene` | Minimum cells expressing a gene |
| `--top_gene_cap` | Number of highly variable genes retained |
| `--umap_init` | UMAP initialization strategy |
| `--umap_min_dist` | UMAP minimum-distance parameter |

---

# Example Outputs

### Marker Scores

```text
cellpose_pipeline_outputs/data/
└── cluster_markers_all_cells_final.csv
```

### Spatial Cell Coordinates

```text
cellpose_pipeline_outputs/data/
└── spatial_clusters_final.csv
```

### Spatial Marker Visualization

```text
cellpose_pipeline_outputs/data/
└── spatial_clusters.csv
```

### UMAP Embeddings

```text
cellpose_pipeline_outputs/data/
├── umap_clusters_final.csv
└── umap_clusters.csv
```

### Seurat Object

```text
cellpose_pipeline_outputs/data/
└── seurat/
```

### Spatial Network

```text
cellpose_pipeline_outputs/data/
├── spatial_graph.edgelist
└── spatial_neighbors.csv
```

### Bayesian Marker Results

```text
cellpose_pipeline_outputs/data/
├── spatial_markers.csv
├── marker_posteriors.csv
└── cluster_markers.csv
```

---

# Benchmarking Framework

SegMap includes benchmarking workflows for comparison against leading spatial transcriptomics analysis pipelines:

- SegMap
- Scanpy
- Seurat
- Baysor

Benchmarking evaluates:

- Spatial coherence
- Moran's I
- Marker specificity
- Marker recovery
- Neighborhood purity
- Silhouette score
- Davies–Bouldin index
- Cluster stability
- Rank robustness

---

# Reproducibility

All analyses presented in the manuscript are fully reproducible.

The repository contains:

- processed datasets
- Synthetic data generation pipelines
- Benchmarking workflows
- Bayesian marker model implementation
- Figure-generation scripts
- R and Python analysis code
- Manuscript source data

---

# Contact

**Dr. Komlan Atitey**  
Department of Neurosciences  
Case Western Reserve University School of Medicine  
komlan.atitey@case.edu
Cleveland, OH, USA

Please use GitHub Issues for bug reports, feature requests, and technical support.
