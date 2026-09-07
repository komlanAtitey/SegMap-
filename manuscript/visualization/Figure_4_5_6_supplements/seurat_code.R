
setwd("/Users/kxa643/Documents/xenium/segmap_R_codes")
getwd()


library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(leidenbase)
# ----------------------------------------------
# 1. LOAD XENIUM DATA (Already a Seurat object)
# ----------------------------------------------

xenium_dir <- "XeniumData/Xenium_V1_FFPE_Human_Brain_Alzheimers_With_Addon_outs"
xenium_dir <- "XeniumData/Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon_outs"
xenium_dir <- "XeniumData/Xenium_V1_FFPE_Human_Brain_Healthy_With_Addon_outs"
xenium_dir <- "XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs"

# Load Xenium data
xenium <- LoadXenium(xenium_dir)

# Check object
Assays(xenium)

# =====================================================
# QC METRICS
# =====================================================

DefaultAssay(xenium) <- "Xenium"

# Total transcripts
xenium$nCount_Xenium <- Matrix::colSums(
  GetAssayData(xenium, layer = "counts")
)

# Detected genes
xenium$nFeature_Xenium <- Matrix::colSums(
  GetAssayData(xenium, layer = "counts") > 0
)

# Mitochondrial percentage (if available)
mt.genes <- grep("^MT-", rownames(xenium), value = TRUE)

if (length(mt.genes) > 0) {
  xenium$percent.mt <- PercentageFeatureSet(
    xenium,
    pattern = "^MT-"
  )
} else {
  xenium$percent.mt <- 0
}

# =====================================================
# VISUALIZE QC
# =====================================================

qc_data <- xenium@meta.data[, c(
  "nFeature_Xenium",
  "nCount_Xenium",
  "percent.mt"
)]

qc_long <- pivot_longer(
  as.data.frame(qc_data),
  cols = everything(),
  names_to = "Metric",
  values_to = "Value"
)

ggplot(
  qc_long,
  aes(Metric, Value, fill = Metric)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    outlier.size = 0.3
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Xenium Quality Control Metrics",
    x = "",
    y = "Value"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

FeatureScatter(
  xenium,
  feature1 = "nCount_Xenium",
  feature2 = "nFeature_Xenium"
)

# =====================================================
# ADAPTIVE QC FILTERING
# (LESS AGGRESSIVE FOR XENIUM)
# =====================================================

summary(xenium$nFeature_Xenium)

feature.low <- max(
  10,
  quantile(
    xenium$nFeature_Xenium,
    0.01,
    na.rm = TRUE
  )
)

feature.high <- quantile(
  xenium$nFeature_Xenium,
  0.995,
  na.rm = TRUE
)

count.low <- max(
  20,
  quantile(
    xenium$nCount_Xenium,
    0.01,
    na.rm = TRUE
  )
)

cat("Feature low cutoff =", feature.low, "\n")
cat("Feature high cutoff =", feature.high, "\n")
cat("Count low cutoff =", count.low, "\n")

xenium <- subset(
  xenium,
  subset =
    nFeature_Xenium > feature.low &
    nFeature_Xenium < feature.high &
    nCount_Xenium > count.low &
    percent.mt < 25
)

cat(
  "Cells retained:",
  ncol(xenium),
  "\n"
)

# =====================================================
# SCTransform
# BETTER FOR XENIUM THAN LOGNORMALIZE
# =====================================================

xenium <- SCTransform(
  xenium,
  assay = "Xenium",
  verbose = FALSE,
  variable.features.n = 3000
)

# =====================================================
# PCA
# =====================================================

xenium <- RunPCA(
  xenium,
  npcs = 50,
  verbose = FALSE
)

ElbowPlot(
  xenium,
  ndims = 50
)

# =====================================================
# NEIGHBOR GRAPH
# =====================================================

dims_use <- 1:30

xenium <- FindNeighbors(
  xenium,
  dims = dims_use,
  k.param = 15
)

# =====================================================
# CLUSTERING
# =====================================================

xenium <- FindClusters(
  xenium,
  resolution = 0.5,
  algorithm = 4
)

# =====================================================
# UMAP
# =====================================================

xenium <- RunUMAP(
  xenium,
  dims = dims_use,
  n.neighbors = 30,
  min.dist = 0.2,
  spread = 1.5,
  metric = "cosine"
)

DimPlot(
  xenium,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

# =====================================================
# 10. SPATIAL VISUALIZATION
# =====================================================
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(xenium)
coords$cluster <- Idents(xenium)

# Custom spatial plot
ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5) +
  scale_y_reverse() +  # y-axis inversion to match Seurat convention
  theme_minimal() +
  labs(title = "Spatial Clusters")

# =====================================================
# MARKERS
# =====================================================
markers <- FindAllMarkers(
  xenium,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.15
)

head(markers)

seurat_data4 <- xenium
save(seurat_data4, file = "seurat_data4.rdata")
seurat_marker4 <- markers
save(seurat_marker4, file = "seurat_marker4.rdata")

# ------------------------------------------------------------------------------
#
#
#                               VISUALIZATION UMAP
#
#
# ------------------------------------------------------------------------------
load("XeniumData_seurat/new_xenium_data1.rdata")
load("XeniumData_seurat/new_xenium_data2.rdata")
load("XeniumData_seurat/new_xenium_data3.rdata")
load("XeniumData_seurat/new_xenium_data4.rdata")

# -------------------------------------------------
#                      UMAP 1
# -------------------------------------------------
DimPlot(
  new_xenium_data1,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "1"="Mature myelinating oligodendrocytes",
  "2"="SPP1⁺ reactive oligodendrocytes",
  "3"="Excitatory glutamatergic neurons (RORB+)",
  "4"="Activated endothelial cells",
  "5"="Astrocytes",
  "6"="Activated microglia / Disease-associated microglia (DAM)",
  "7"="GABAergic interneurons",
  "8"="Oligodendrocyte precursor cells (OPCs / NG2 glia)"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
new_xenium_data1 <- RenameIdents(
  new_xenium_data1,
  cluster_labels
)

# UMAP
p_seurat <- DimPlot(
  new_xenium_data1,
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
p_seurat
p_seurat + NoLegend()

# -------------------------------------------------
#                      SPATIAL 1
# -------------------------------------------------
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(new_xenium_data1)
coords$cluster <- Idents(new_xenium_data1)

# Convert cluster to factor with desired order
coords$cluster <- factor(
  coords$cluster,
  levels = names(cluster_labels),
  labels = cluster_labels
)

p_seurat_spatial <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = setNames(cluster_colors, cluster_labels),
    name = "Cell type"
  ) +
  scale_y_reverse() +   # Match Seurat/Xenium orientation
  coord_fixed() +       # Preserve spatial aspect ratio
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_seurat_spatial
p_seurat_spatial + NoLegend()

# --------------------------------
#               UMAP 2
# --------------------------------
DimPlot(
  new_xenium_data2,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "1"="Mature myelinating oligodendrocytes",
  "2"="Vascular endothelial cells",
  "3"="Activated microglia / disease-associated microglia (DAM)",
  "4"="HTR2A⁺ excitatory glutamatergic neurons",
  "5"="Excitatory glutamatergic neurons",
  "6"="Reactive astrocytes",
  "7"="Oligodendrocyte precursor cells (OPCs / NG2 glia)",
  "8"="LHX6⁺ inhibitory interneurons",
  "9"="Excitatory projection neurons",
  "10"="CD8⁺ T lymphocytes",
  "11"="Deep-layer corticothalamic excitatory neurons"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
new_xenium_data2 <- RenameIdents(
  new_xenium_data2,
  cluster_labels
)

# UMAP
p_seurat <- DimPlot(
  new_xenium_data2,
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
p_seurat
p_seurat + NoLegend()


# -------------------------------------------------
#                      SPATIAL 2
# -------------------------------------------------
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(new_xenium_data2)
coords$cluster <- Idents(new_xenium_data2)

# Convert cluster to factor with desired order
coords$cluster <- factor(
  coords$cluster,
  levels = names(cluster_labels),
  labels = cluster_labels
)

p_seurat_spatial <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = setNames(cluster_colors, cluster_labels),
    name = "Cell type"
  ) +
  scale_y_reverse() +   # Match Seurat/Xenium orientation
  coord_fixed() +       # Preserve spatial aspect ratio
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_seurat_spatial
p_seurat_spatial + NoLegend()

# --------------------------------
#               UMAP 3
# --------------------------------
DimPlot(
  new_xenium_data3,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "1"="LAMP5⁺ upper-layer excitatory neurons (CUX2⁺ intratelencephalic neurons)",
  "2"="Mature myelinating oligodendrocytes",
  "3"="Reactive astrocytes",
  "4"="Vascular endothelial cells",
  "5"="LHX6⁺ inhibitory interneurons (PV/SST lineage)",
  "6"="Excitatory glutamatergic projection neurons",
  "7"="Activated microglia / Disease-associated microglia (DAM)",
  "8"="Perivascular fibroblasts / vascular fibroblasts",
  "9"="VIP⁺ inhibitory interneurons",
  "10"="Oligodendrocyte precursor cells (OPCs / NG2 glia)",
  "11"="RELN⁺ inhibitory interneurons (neurogliaform/LAMP5-like)"
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
new_xenium_data3 <- RenameIdents(
  new_xenium_data3,
  cluster_labels
)

# UMAP
p_seurat <- DimPlot(
  new_xenium_data3,
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
p_seurat
p_seurat + NoLegend()

# -------------------------------------------------
#                      SPATIAL 3
# -------------------------------------------------
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(new_xenium_data3)
coords$cluster <- Idents(new_xenium_data3)

# Convert cluster to factor with desired order
coords$cluster <- factor(
  coords$cluster,
  levels = names(cluster_labels),
  labels = cluster_labels
)

p_seurat_spatial <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = setNames(cluster_colors, cluster_labels),
    name = "Cell type"
  ) +
  scale_y_reverse() +   # Match Seurat/Xenium orientation
  coord_fixed() +       # Preserve spatial aspect ratio
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_seurat_spatial
p_seurat_spatial + NoLegend()

# --------------------------------
#               UMAP 4
# --------------------------------
DimPlot(
  new_xenium_data4,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Xenium Cell Clusters")

#===========
cluster_labels <- c(
  "1"="Activated / memory CD4⁺ T cells",
  "2"="Epithelial cells (secretory/luminal epithelial)",
  "3"="Blood vascular endothelial cells",
  "4"="Activated macrophages (M2-like tissue macrophages)",
  "5"="Matrix-producing fibroblasts",
  "6"="Vascular smooth muscle cells",
  "7"="Ciliated epithelial cells",
  "8"="B cells",
  "9"="Alveolar type II (AT2) epithelial cells",
  "10"="Plasma cells",
  "11"="Alveolar type I (AT1) epithelial cells",
  "12"="Mast cells",
  "13"="Lymphatic endothelial cells"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")


# Rename cluster identities
new_xenium_data4 <- RenameIdents(
  new_xenium_data4,
  cluster_labels
)

# UMAP
p_seurat <- DimPlot(
  new_xenium_data4,
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
p_seurat
p_seurat + NoLegend()

# -------------------------------------------------
#                      SPATIAL 4
# -------------------------------------------------
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(new_xenium_data4)
coords$cluster <- Idents(new_xenium_data4)

# Convert cluster to factor with desired order
coords$cluster <- factor(
  coords$cluster,
  levels = names(cluster_labels),
  labels = cluster_labels
)

p_seurat_spatial <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = setNames(cluster_colors, cluster_labels),
    name = "Cell type"
  ) +
  scale_y_reverse() +   # Match Seurat/Xenium orientation
  coord_fixed() +       # Preserve spatial aspect ratio
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_seurat_spatial
p_seurat_spatial + NoLegend()











