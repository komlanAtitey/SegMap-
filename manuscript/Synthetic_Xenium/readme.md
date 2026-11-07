

# Synthetic Xenium Data Generation and SegMap Benchmarking

## Overview

This repository provides a workflow for generating synthetic spatial transcriptomics data, converting the synthetic tissue into a Xenium-compatible dataset, running the SegMap segmentation pipeline, and evaluating segmentation performance against the synthetic ground truth.

The workflow consists of five main steps:

1. Generate synthetic tissue containing 900 cells per cell type.
2. Render the synthetic tissue into a Xenium-compatible directory.
3. Run SegMap segmentation on the synthetic Xenium dataset.
4. Convert the SegMap segmentation mask into a prediction CSV file.
5. Evaluate segmentation performance against the synthetic ground truth.

The primary objective is to establish a reproducible benchmarking framework for assessing SegMap segmentation performance using synthetic spatial transcriptomics data with known ground-truth annotations.

---

## 1. Project Structure

The expected project structure is:

```text
project_directory/
│
├── segmap_synthetic_GitHub/
│   ├── synthetic_tissue.py
│   ├── synth_to_xenium.py
│   ├── segmap_to_pred.py
│   └── score_pipeline.py
│
├── synthetic_out/
│
├── synth_xenium/
│   └── segmap_out/
│
├── segmap_pred.csv
│
└── score_segmap/
```

The output directories and files are generated during execution.

---

## 2. Requirements

Before running the workflow, ensure that the following requirements are met:

- Python is installed and accessible from the terminal.
- SegMap is installed in the active Python environment.
- Required Python dependencies for the synthetic data generation scripts are installed.
- GPU support is configured if using the `--gpu` option.
- Sufficient disk space is available for synthetic tissue generation, Xenium rendering, segmentation masks, and evaluation outputs.

Activate the appropriate Python environment before executing the commands.

Verify that SegMap is accessible:

```bash
python -m segmap.run_segmap --help
```

---

## 3. Step 1: Generate Synthetic Tissue

Generate synthetic tissue containing 900 cells per cell type.

```bash
python segmap_synthetic_GitHub/synthetic_tissue.py \
    --n-cells 900 \
    --seed 0 \
    --out-dir synthetic_out
```

### Parameters

| Parameter | Description | Value |
|---|---|---|
| `--n-cells` | Number of cells per cell type | 900 |
| `--seed` | Random seed for reproducibility | 0 |
| `--out-dir` | Output directory for synthetic ground-truth data | `synthetic_out` |

### Output

The generated synthetic tissue and associated ground-truth data are saved in:

```text
synthetic_out/
```

This directory serves as the reference dataset for subsequent segmentation evaluation.

---

## 4. Step 2: Render Synthetic Tissue into Xenium Format

Convert the generated synthetic tissue into a Xenium-compatible directory.

This step renders the synthetic tissue and generates a new DAPI image corresponding to the synthetic nuclei.

```bash
python segmap_synthetic_GitHub/synth_to_xenium.py \
    --gt-dir synthetic_out \
    --out-dir synth_xenium \
    --scale 2
```

### Parameters

| Parameter | Description | Value |
|---|---|---|
| `--gt-dir` | Synthetic ground-truth directory | `synthetic_out` |
| `--out-dir` | Xenium-compatible output directory | `synth_xenium` |
| `--scale` | Rendering scale factor | 2 |

### Output

The Xenium-compatible synthetic dataset is saved in:

```text
synth_xenium/
```

This directory is used as the input for SegMap.

---

## 5. Step 3: Run SegMap Segmentation

Run SegMap on the newly rendered synthetic Xenium dataset.

### 5.1. Remove Previous SegMap Results

Remove the previous SegMap output directory before starting a new run:

```bash
rm -rf synth_xenium/segmap_out
```

**Important:** This command permanently removes the existing `segmap_out` directory and its contents. Execute it only when you intend to discard the previous segmentation results.

This ensures that the new run does not reuse the previous segmentation mask, including the previously generated 2477 × 2477 mask.

### 5.2. Execute SegMap

Run the segmentation pipeline:

```bash
python -m segmap.run_segmap \
    --xenium_dir synth_xenium \
    --outdir synth_xenium/segmap_out \
    --no_tune_resolution \
    --resolution 0.6 \
    --n_pcs 20 \
    --cellpose_diameter 23 \
    --expand_pixels 15 \
    --min_genes_per_cell 10 \
    --min_counts_per_cell 25 \
    --min_cells_per_gene 3 \
    --gpu
```

### SegMap Parameters

| Parameter | Description | Value |
|---|---|---|
| `--xenium_dir` | Input Xenium dataset | `synth_xenium` |
| `--outdir` | SegMap output directory | `synth_xenium/segmap_out` |
| `--no_tune_resolution` | Disable automatic resolution tuning | Enabled |
| `--resolution` | Clustering resolution | 0.6 |
| `--n_pcs` | Number of principal components | 20 |
| `--cellpose_diameter` | Cellpose diameter parameter | 23 |
| `--expand_pixels` | Cell boundary expansion parameter | 15 |
| `--min_genes_per_cell` | Minimum genes required per cell | 10 |
| `--min_counts_per_cell` | Minimum transcript counts required per cell | 25 |
| `--min_cells_per_gene` | Minimum cells expressing a gene | 3 |
| `--gpu` | Enable GPU execution | Enabled |

### Output

SegMap results are saved in:

```text
synth_xenium/segmap_out/
```

The segmentation mask used in the subsequent evaluation is:

```text
synth_xenium/segmap_out/data/full_mask_cellpose_fullres_tiled.tif
```

---

## 6. Step 4: Convert SegMap Mask into Prediction CSV

Convert the SegMap segmentation mask into a prediction CSV file for comparison with the synthetic ground truth.

```bash
python segmap_synthetic_GitHub/segmap_to_pred.py \
    --gt-dir synthetic_out \
    --mask synth_xenium/segmap_out/data/full_mask_cellpose_fullres_tiled.tif \
    --out segmap_pred.csv
```

### Parameters

| Parameter | Description |
|---|---|
| `--gt-dir` | Synthetic ground-truth directory |
| `--mask` | SegMap-generated segmentation mask |
| `--out` | Output prediction CSV file |

### Output

The predicted segmentation results are saved as:

```text
segmap_pred.csv
```

This file contains the SegMap-derived predictions used for downstream evaluation.

---

## 7. Step 5: Evaluate Segmentation Performance

Evaluate the SegMap predictions against the synthetic ground-truth dataset.

```bash
python segmap_synthetic_GitHub/score_pipeline.py \
    --gt-dir synthetic_out \
    --pred segmap_pred.csv \
    --out-dir score_segmap
```

### Parameters

| Parameter | Description |
|---|---|
| `--gt-dir` | Synthetic ground-truth directory |
| `--pred` | SegMap prediction CSV file |
| `--out-dir` | Directory for evaluation results |

### Output

Evaluation results are saved in:

```text
score_segmap/
```

These results provide the basis for assessing SegMap segmentation performance against the known synthetic ground truth.

The specific evaluation metrics and output files depend on the implementation of `score_pipeline.py`.

---

## 8. Complete Workflow: Run All Steps

The entire workflow can also be executed sequentially using a single shell script.

Create a file named:

```text
run_synthetic_segmap.sh
```

Add the following commands:

```bash
#!/bin/bash

set -e

echo "=========================================="
echo "STEP 1: Generate synthetic tissue"
echo "=========================================="

python segmap_synthetic_GitHub/synthetic_tissue.py \
    --n-cells 900 \
    --seed 0 \
    --out-dir synthetic_out


echo "=========================================="
echo "STEP 2: Render Xenium-compatible dataset"
echo "=========================================="

python segmap_synthetic_GitHub/synth_to_xenium.py \
    --gt-dir synthetic_out \
    --out-dir synth_xenium \
    --scale 2


echo "=========================================="
echo "STEP 3: Run SegMap"
echo "=========================================="

rm -rf synth_xenium/segmap_out

python -m segmap.run_segmap \
    --xenium_dir synth_xenium \
    --outdir synth_xenium/segmap_out \
    --no_tune_resolution \
    --resolution 0.6 \
    --n_pcs 20 \
    --cellpose_diameter 23 \
    --expand_pixels 15 \
    --min_genes_per_cell 10 \
    --min_counts_per_cell 25 \
    --min_cells_per_gene 3 \
    --gpu


echo "=========================================="
echo "STEP 4: Convert SegMap mask to predictions"
echo "=========================================="

python segmap_synthetic_GitHub/segmap_to_pred.py \
    --gt-dir synthetic_out \
    --mask synth_xenium/segmap_out/data/full_mask_cellpose_fullres_tiled.tif \
    --out segmap_pred.csv


echo "=========================================="
echo "STEP 5: Evaluate segmentation performance"
echo "=========================================="

python segmap_synthetic_GitHub/score_pipeline.py \
    --gt-dir synthetic_out \
    --pred segmap_pred.csv \
    --out-dir score_segmap


echo "=========================================="
echo "WORKFLOW COMPLETED SUCCESSFULLY"
echo "=========================================="
```

Make the script executable:

```bash
chmod +x run_synthetic_segmap.sh
```

Execute the complete workflow:

```bash
./run_synthetic_segmap.sh
```

---

## 9. Expected Final Outputs

After successful execution, the main outputs should be organized as follows:

| Output | Description |
|---|---|
| `synthetic_out/` | Synthetic tissue and ground-truth data |
| `synth_xenium/` | Xenium-compatible synthetic dataset |
| `synth_xenium/segmap_out/` | SegMap segmentation results |
| `full_mask_cellpose_fullres_tiled.tif` | SegMap segmentation mask |
| `segmap_pred.csv` | SegMap-derived prediction data |
| `score_segmap/` | Segmentation evaluation results |

---

## 10. Reproducibility and Benchmarking Considerations

The synthetic data generation uses a fixed random seed (`--seed 0`) to support reproducibility.

For consistent benchmarking:

1. Keep the synthetic tissue generation parameters fixed across experiments.
2. Use the same synthetic ground-truth dataset when comparing segmentation methods.
3. Maintain consistent Xenium rendering parameters.
4. Record SegMap configuration parameters for each experiment.
5. Preserve the generated segmentation masks and prediction CSV files.
6. Store evaluation results separately for each segmentation method or parameter configuration.

For additional experiments, the number of cells, random seed, rendering scale, and SegMap parameters can be modified systematically.

---

## 11. Workflow Summary

```text
Synthetic Tissue Generation
          |
          v
     synthetic_out/
          |
          v
  Xenium Data Rendering
          |
          v
     synth_xenium/
          |
          v
     SegMap Pipeline
          |
          v
    Segmentation Mask
          |
          v
    Mask-to-Prediction
          |
          v
     segmap_pred.csv
          |
          v
   Ground-Truth Scoring
          |
          v
      score_segmap/
```

**Final objective:** Generate reproducible synthetic Xenium datasets, execute SegMap segmentation, and quantitatively evaluate its segmentation predictions against known synthetic ground-truth annotations.
