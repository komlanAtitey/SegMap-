setwd("/Users/kxa643/Documents/xenium/segmap_R_codes")
getwd()

###############################################################
## Baysor analysis pipeline for 10x Xenium
## Author: Komlan Atitey
###############################################################

rm(list=ls())

library(Seurat)
library(data.table)
library(Matrix)
library(dplyr)
library(R.utils)
library(ggplot2)

###############################################################
## USER SETTINGS
###############################################################

## Xenium output directory
xenium.dir <- "XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs"

## Baysor executable
baysor <- "Baysor"

## Output directory
outdir <- file.path(xenium.dir,"Baysor")

dir.create(outdir,showWarnings=FALSE)

###############################################################
## Read transcript table
###############################################################

cat("Reading transcripts...\n")

transcripts <- fread(
  file.path(
    xenium.dir,
    "transcripts.csv.gz"
  )
)

#-------------------------------------------------------
# In the case there is no .csv.gz, use the block below:
#-------------------------------------------------------
library(arrow)
transcripts <- read_parquet(
  file.path(
    xenium.dir,
    "transcripts.parquet"
  )
)

transcripts <- as.data.table(transcripts)

###############################################################
## Filter transcripts
###############################################################

cat("Filtering Q20 transcripts...\n")

transcripts <- transcripts[
  qv >= 20 &
    feature_name != "NegControlProbe" &
    feature_name != "NegControlCodeword"
]

###############################################################
## Remove unassigned transcripts
###############################################################

transcripts <- transcripts[
  !is.na(cell_id)
]

###############################################################
## Save Baysor input
###############################################################

input.file <- file.path(
  outdir,
  "baysor_input.csv"
)

fwrite(
  transcripts,
  input.file
)

###############################################################
## Run Baysor
###############################################################

cat("Running Baysor...\n")

cmd <- paste(
  
  baysor,
  
  "run",
  
  "-x x_location",
  "-y y_location",
  "-z z_location",
  
  "-g feature_name",
  
  "-m 20",
  
  "--scale 15",
  
  "--prior-segmentation-confidence 0.3",
  
  "-o", outdir,
  
  input.file,
  
  ":cell_id"
  
)

system(cmd)

#========
baysor_exec <- "/Users/kxa643/Downloads/Baysor/build/user/baysor"
experiment <- "/Users/kxa643/Documents/xenium/XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs/experiment.xenium"
config <- "/Users/kxa643/Downloads/Baysor/configs/xenium.toml"
outdir <- "/Users/kxa643/Documents/xenium/XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs/Baysor"

res <- system2(
  baysor_exec,
  args = c(
    "run",
    experiment,
    ":cell_id",
    "-c", config,
    "-o", outdir
  ),
  stdout = TRUE,
  stderr = TRUE
)

cat(res, sep = "\n")
status <- attr(res, "status")

if (!is.null(status)) {
  stop("Baysor failed.")
}

###############################################################
## Check Baysor output
###############################################################

print(list.files(outdir))

seg.file <- file.path(outdir, "segmentation.csv")

if (!file.exists(seg.file)) {
  
  seg.file <- list.files(
    outdir,
    pattern = "segment",
    full.names = TRUE
  )[1]
  
  if (is.na(seg.file))
    stop("No segmentation file produced by Baysor.")
}

cat("Using segmentation file:\n")
cat(seg.file, "\n")

###############################################################
## Read Baysor segmentation
###############################################################
library(data.table)

seg <- fread(seg.file)

###############################################################
## Detect cell column
###############################################################

if ("cell" %in% names(seg)) {
  
  cell.col <- "cell"
  
} else if ("cell_id" %in% names(seg)) {
  
  cell.col <- "cell_id"
  
} else {
  
  stop("Cannot identify cell column.")
}

###############################################################
## Detect gene column
###############################################################

if ("gene" %in% names(seg)) {
  
  gene.col <- "gene"
  
} else if ("feature_name" %in% names(seg)) {
  
  gene.col <- "feature_name"
  
} else {
  
  stop("Cannot identify gene column.")
}

###############################################################
## Remove noise
###############################################################

if ("is_noise" %in% names(seg)) {
  seg <- seg[is_noise == FALSE]
  
}

###############################################################
## Build count matrix
###############################################################

counts <- seg[
  ,
  .N,
  by = c(gene.col, cell.col)
]

setnames(
  counts,
  c(gene.col, cell.col),
  c("gene", "cell")
)

counts.mat <- dcast(
  counts,
  gene ~ cell,
  value.var = "N",
  fill = 0
)

genes <- counts.mat$gene

counts.mat <- as.matrix(
  counts.mat[, -1]
)

rownames(counts.mat) <- genes

counts.mat <- Matrix::Matrix(
  counts.mat,
  sparse = TRUE
)

###############################################################
## Compute centroids
###############################################################
coord.cols <- intersect(
  c("x", "y", "z"),
  names(seg)
)

coords <- seg[
  ,
  lapply(.SD, mean),
  by = cell.col,
  .SDcols = coord.cols
]

setnames(coords, cell.col, "cell")
coords <- as.data.frame(coords)
rownames(coords) <- coords$cell

###############################################################
## Create Seurat object
###############################################################

obj <- CreateSeuratObject(
  counts=counts.mat,
  assay="RNA",
  project="Baysor"
)

###############################################################
## Compute cell centroids
###############################################################

coords <- seg[
  ,
  .(
    x=mean(x),
    y=mean(y),
    z=mean(z)
  ),
  by=cell
]

coords <- as.data.frame(coords)

rownames(coords) <- coords$cell

obj <- AddMetaData(
  obj,
  metadata=coords
)

###############################################################
## Standard Seurat workflow
###############################################################

obj <- NormalizeData(obj)

obj <- FindVariableFeatures(obj)

obj <- ScaleData(obj)

obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims=1:30)
obj <- FindClusters(obj, resolution=0.16)
obj <- RunUMAP(obj, dims=1:30)

DimPlot(
  obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) 

###############################################################
# Marker genes
###############################################################
markers <- FindAllMarkers(
  obj,
  only.pos=TRUE,
  min.pct=.25,
  logfc.threshold=.25
)

baysor_marker4 <- markers
baysor_data_4 <- obj


# =================================================
save(baysor_marker4, file="baysor_marker4.rdata")
save(baysor_data_4, file="baysor_data_4.rdata")


# ------------------------------------------------------------------------------
#
#
#                               VISUALIZATION UMAP
#
#
# ------------------------------------------------------------------------------
load("XeniumData_baysor/baysor_marker1.rdata")
load("XeniumData_baysor/baysor_marker2.rdata")
load("XeniumData_baysor/baysor_marker3.rdata")
load("XeniumData_baysor/baysor_marker4.rdata")

marker_1 <- baysor_marker1[,6:7]
marker_2 <- baysor_marker2[,6:7]
marker_3 <- baysor_marker3[,6:7]
marker_4 <- baysor_marker4[,6:7]

write.csv(
  marker_4,
  file = "marker_4.csv",
  row.names = FALSE
)

load("XeniumData_baysor/baysor_data_1.rdata")
load("XeniumData_baysor/baysor_data_2.rdata")
load("XeniumData_baysor/baysor_data_3.rdata")
load("XeniumData_baysor/baysor_data_4.rdata")

# -------------------------------------------------
#                      UMAP 1
# -------------------------------------------------
DimPlot(
  baysor_data_1,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "0"="Homeostatic astrocytes",
  "1"="Reactive oligodendrocytes",
  "2"="Excitatory neurons",
  "3"="Activated endothelial cells",
  "4"="Mixed (astrocyte + oligodendrocyte + neuron)",
  "5"="OPCs (oligodendrocyte precursor cells)",
  "6"="Stressed neurons",
  "7"="Apoptotic / severely stressed cells"
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
baysor_data_1 <- RenameIdents(
  baysor_data_1,
  cluster_labels
)

# UMAP
p_baysor <- DimPlot(
  baysor_data_1,
  reduction = "umap",
  cols = cluster_colors,
  label = FALSE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.5
) +
  ggtitle("Annotated Xenium Cell Types") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    )
  )
p_baysor
p_baysor + NoLegend()


# -------------------------------------------------
#                      SPATIAL 1
# -------------------------------------------------
celltype_colors <- setNames(
  cluster_colors,
  unname(cluster_labels)
)

## Metadata
spatial.df <- baysor_data_1@meta.data
spatial.df$cell <- rownames(spatial.df)
spatial.df$cluster <- as.character(spatial.df$seurat_clusters)

## Add cell-type labels
spatial.df$cell_type <- cluster_labels[spatial.df$cluster]

## Remove cells without annotations (optional)
spatial.df <- spatial.df[!is.na(spatial.df$cell_type), ]

## Keep legend in desired order
spatial.df$cell_type <- factor(
  spatial.df$cell_type,
  levels = unname(cluster_labels)
)

## Spatial plot
p_baysor_spatial <- ggplot(
  spatial.df,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(size = 0.35, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    title = "Baysor Spatial Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_baysor_spatial
p_baysor_spatial + NoLegend()

# --------------------------------
#               UMAP 2
# --------------------------------
DimPlot(
  baysor_data_2,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "0"="L2/3–L4 Intratelencephalic (IT) excitatory neurons (RORB⁺)",
  "1"="Astrocytes (reactive / fibrous)",
  "2"="Mature myelinating oligodendrocytes",
  "3"="Activated microglia / disease-associated microglia (DAM)",
  "4"="Endothelial cells / vascular endothelial cells",
  "5"="OPCs (oligodendrocyte precursor cells)",
  "6"="SST/PVALB interneurons (LHX6⁺ MGE-derived inhibitory neurons)",
  "7"="Deep-layer corticocortical excitatory neurons (L6 IT / corticothalamic-like)",
  "8"="VIP/LAMP5 interneurons",
  "9"="T lymphocytes (with some macrophage contamination)",
  "10"="LAMP5 neurogliaform interneurons)"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
baysor_data_2 <- RenameIdents(
  baysor_data_2,
  cluster_labels
)

# UMAP
p_baysor <- DimPlot(
  baysor_data_2,
  reduction = "umap",
  cols = cluster_colors,
  label = FALSE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.5
) +
  ggtitle("Annotated Xenium Cell Types") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    )
  )

p_baysor
p_baysor + NoLegend()


# -------------------------------------------------
#                      SPATIAL 2
# -------------------------------------------------
celltype_colors <- setNames(
  cluster_colors,
  unname(cluster_labels)
)

## Metadata
spatial.df <- baysor_data_2@meta.data
spatial.df$cell <- rownames(spatial.df)
spatial.df$cluster <- as.character(spatial.df$seurat_clusters)

## Add cell-type labels
spatial.df$cell_type <- cluster_labels[spatial.df$cluster]

## Remove cells without annotations (optional)
spatial.df <- spatial.df[!is.na(spatial.df$cell_type), ]

## Keep legend in desired order
spatial.df$cell_type <- factor(
  spatial.df$cell_type,
  levels = unname(cluster_labels)
)

## Spatial plot
p_baysor_spatial <- ggplot(
  spatial.df,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(size = 0.35, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    title = "Baysor Spatial Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_baysor_spatial
p_baysor_spatial + NoLegend()


# --------------------------------
#               UMAP 3
# --------------------------------
DimPlot(
  baysor_data_3,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "0"="Excitatory neurons (L2/3 IT, RORB+)",
  "1"="Excitatory neurons (L2/3 IT, LAMP5/CCK subtype)",
  "2"="Mature oligodendrocytes",
  "3"="Astrocytes",
  "4"="Endothelial cells",
  "5"="OPVALB interneurons",
  "6"="VIP/LAMP5 interneurons",
  "7"="Microglia",
  "8"="Neuronal doublets / mixed cells (review recommended)",
  "9"="Oligodendrocyte precursor cells (OPCs)",
  "10"="Excitatory neurons (L2/3 IT, NPNT/MCTP2 subtype)"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# Rename cluster identities
baysor_data_3 <- RenameIdents(
  baysor_data_3,
  cluster_labels
)

# UMAP
p_baysor <- DimPlot(
  baysor_data_3,
  reduction = "umap",
  cols = cluster_colors,
  label = FALSE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.5
) +
  ggtitle("Annotated Xenium Cell Types") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    )
  )
p_baysor
p_baysor + NoLegend()

# -------------------------------------------------
#                      SPATIAL 3
# -------------------------------------------------
celltype_colors <- setNames(
  cluster_colors,
  unname(cluster_labels)
)

## Metadata
spatial.df <- baysor_data_3@meta.data
spatial.df$cell <- rownames(spatial.df)
spatial.df$cluster <- as.character(spatial.df$seurat_clusters)

## Add cell-type labels
spatial.df$cell_type <- cluster_labels[spatial.df$cluster]

## Remove cells without annotations (optional)
spatial.df <- spatial.df[!is.na(spatial.df$cell_type), ]

## Keep legend in desired order
spatial.df$cell_type <- factor(
  spatial.df$cell_type,
  levels = unname(cluster_labels)
)

## Spatial plot
p_baysor_spatial <- ggplot(
  spatial.df,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(size = 0.35, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    title = "Baysor Spatial Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_baysor_spatial
p_baysor_spatial + NoLegend()

# --------------------------------
#               UMAP 4
# --------------------------------
DimPlot(
  baysor_data_4,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "0"="Malignant epithelial cells (LUAD, EPCAM+/EGFR+)",
  "1"="Tumor-associated macrophages (TAMs, M2-like)",
  "2"="T lymphocytes (activated/memory CD4/CD8)",
  "3"="Cancer-associated fibroblasts (myofibroblastic CAFs)",
  "4"="Alveolar capillary endothelial cells",
  "5"="Venous endothelial cells",
  "6"="Alveolar type II epithelial cells (AT2) / LUAD epithelial cells",
  "7"="Pericytes / vascular smooth muscle cells",
  "8"="B cells / plasma cells",
  "9"="Secretory/ciliated tumor epithelial cells (club/goblet-like)",
  "10"="Proliferating malignant epithelial cells",
  "11"="Mast cells",
  "12"="Lymphatic endothelial cells"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# Rename cluster identities
baysor_data_4 <- RenameIdents(
  baysor_data_4,
  cluster_labels
)

# UMAP
p_baysor <- DimPlot(
  baysor_data_4,
  reduction = "umap",
  cols = cluster_colors,
  label = FALSE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.5,
  raster = FALSE
) +
  ggtitle("Annotated Xenium Cell Types") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    )
  )
p_baysor
p_baysor + NoLegend()

# -------------------------------------------------
#                      SPATIAL 4
# -------------------------------------------------
celltype_colors <- setNames(
  cluster_colors,
  unname(cluster_labels)
)

## Metadata
spatial.df <- baysor_data_4@meta.data
spatial.df$cell <- rownames(spatial.df)
spatial.df$cluster <- as.character(spatial.df$seurat_clusters)

## Add cell-type labels
spatial.df$cell_type <- cluster_labels[spatial.df$cluster]

## Remove cells without annotations (optional)
spatial.df <- spatial.df[!is.na(spatial.df$cell_type), ]

## Keep legend in desired order
spatial.df$cell_type <- factor(
  spatial.df$cell_type,
  levels = unname(cluster_labels)
)

## Spatial plot
p_baysor_spatial <- ggplot(
  spatial.df,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(size = 0.35, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  coord_fixed() +
  scale_y_reverse() +
  labs(
    title = "Baysor Spatial Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_baysor_spatial
p_baysor_spatial + NoLegend()
















































