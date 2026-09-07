setwd("/Users/kxa643/Documents/xenium/segmap_R_codes")
getwd()

#===================================================
#               Load required packages
#===================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(leidenbase)
library(patchwork)
library(Matrix)
library(ggplot2)
library(readr)
library(tibble)
library(knitr)
library(purrr)
library(spdep)

library(reticulate)
py_config()
conda_list()

Sys.setenv(
  RETICULATE_PYTHON = "/Users/kxa643/anaconda3/bin/python"
)
use_python(
  "/Users/kxa643/anaconda3/bin/python",
  required = TRUE
)
py_config()
#===================================================
#               Read Xenium matrix with Scanpy
#===================================================
sc <- import("scanpy")

adata <- sc$read_10x_h5(
  "/Users/kxa643/Documents/xenium/XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs/cell_feature_matrix.h5" 
)

#===================================================
#               Scanpy preprocessing
#===================================================
sc$pp$normalize_total(adata, target_sum = 1e4)
sc$pp$log1p(adata)
sc$tl$pca(adata)
sc$pp$neighbors(adata)
sc$tl$leiden(adata, resolution = 0.82, flavor = "igraph", directed = FALSE, n_iterations = 2L)
sc$tl$umap(adata)

#===================================================
#               UMAP visualization
#===================================================
# Display UMAP
sc$pl$umap(
  adata,
  color = "leiden",
  legend_loc = "on data",
  frameon = FALSE,
  size = 8,
  show = TRUE
)

# Save high-resolution UMAP
sc$settings$figdir <- "scanpy_data"

sc$pl$umap(
  adata,
  color = "leiden",
  legend_loc = "on data",
  frameon = FALSE,
  size = 8,
  show = FALSE,
  save = "_leiden_clusters.png"
)

#===================================================
#               Create output directory
#===================================================
dir.create("scanpy_data", showWarnings = FALSE)

#===================================================
#               Export count matrix
#===================================================

counts_mat <- as.matrix(adata$X)
cell_names <- py_to_r(adata$obs_names$to_list())
gene_names <- py_to_r(adata$var_names$to_list())
counts_df <- as.data.frame(counts_mat)
rownames(counts_df) <- cell_names
colnames(counts_df) <- gene_names
write.csv( counts_df,"scanpy_data/cell_by_gene.csv")

#===================================================
#               Export UMAP
#===================================================
# Extract UMAP coordinates
umap_df <- as.data.frame(py_to_r(adata$obsm["X_umap"]))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
rownames(umap_df) <- rownames(adata$obs)

# Add Leiden clusters
umap_df$cluster <- factor(adata$obs$leiden)
write.csv(umap_df, "scanpy_data/umap_with_clusters.csv", row.names = TRUE)

#===================================================
#               Publication-quality UMAP in R
#===================================================

umap_plot <- ggplot(umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.6, alpha = 0.8) +
  theme_classic(base_size = 14) +
  labs(title = "Scanpy Leiden Clusters", x = "UMAP 1", y = "UMAP 2", color = "Cluster") +
  theme(aspect.ratio = 1, legend.position = "right" )
print(umap_plot)

ggsave("scanpy_data/umap_leiden_clusters_R.pdf",
  umap_plot,width = 8, height = 7, dpi = 600)

ggsave("scanpy_data/umap_leiden_clusters_R.png",
  umap_plot,width = 8, height = 7, dpi = 600)

#===================================================
#               Extract Leiden clusters
#===================================================
leiden_clusters <- py_to_r(adata$obs["leiden"])$leiden
umap_df$cluster <- leiden_clusters
write.csv(umap_df,"scanpy_data/umap_with_clusters.csv")


#===================================================
#               Load Xenium object
#===================================================
xenium_dir <-"XeniumData/Xenium_V1_Human_Lung_Cancer_FFPE_outs"
xenium <- LoadXenium(xenium_dir)

#===================================================
#               Extract spatial coordinates
#===================================================
spatial_coords <- GetTissueCoordinates(xenium)
rownames(spatial_coords) <- colnames(xenium)

#---------------------------------------------------
# Ensure same cell order
#---------------------------------------------------
if (!all(rownames(spatial_coords) == cell_names)) {
  spatial_coords <-
    spatial_coords[cell_names, ]
}

write.csv(spatial_coords, "scanpy_data/spatial_coords.csv")
adata$obsm["spatial"] <-r_to_py(as.matrix(spatial_coords))

#===================================================
#               Spatial visualization
#===================================================

spatial_df <- as.data.frame(spatial_coords)
colnames(spatial_df)[1:2] <- c("x", "y")
spatial_df$cluster <- factor(adata$obs$leiden)

spatial_plot <- ggplot(
    spatial_df, aes(x = x, y = y, color = cluster)) +
  geom_point( size = 0.35, alpha = 0.9) +
  scale_y_reverse() +
  coord_fixed() +
  theme_classic(base_size = 14) +
  labs(title = "Spatial Distribution of Scanpy Leiden Clusters",
    x = "X coordinate", y = "Y coordinate",color = "Cluster")
print(spatial_plot)

ggsave("scanpy_data/spatial_clusters.pdf",
  spatial_plot,width = 8,height = 8,dpi = 600)

ggsave("scanpy_data/spatial_clusters.png",
       spatial_plot,width = 8,height = 8,dpi = 600)

#===================================================
#               Transfer Scanpy clusters
#===================================================
if (!all(colnames(xenium) == cell_names)) {
  xenium <- xenium[, cell_names]}

xenium$scanpy_cluster <- factor(leiden_clusters)
Idents(xenium) <- "scanpy_cluster"

#===================================================
#               Marker genes (Scanpy)
#===================================================
# Wilcoxon differential expression
sc$tl$rank_genes_groups(
  adata,
  groupby = "leiden",
  method = "wilcoxon"
)

# Convert results to an R data frame
markers <- py_to_r(
  sc$get$rank_genes_groups_df(
    adata,
    group = NULL
  )
)

write.csv(
  markers,
  "scanpy_data/marker_genes.csv",
  row.names = FALSE
)

message("Marker genes saved.")

write.csv(
  markers,
  "scanpy_data/marker_genes.csv",
  row.names = FALSE
)

message(
  "Marker genes saved."
)

#===================================================
#               Top markers
#===================================================
library(dplyr)

top_markers <-
  markers %>%
  group_by(group) %>%
  slice_max(
    order_by = logfoldchanges,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top_markers,
  "scanpy_data/top_marker_genes_per_cluster.csv",
  row.names = FALSE
)

message("Top marker genes saved.")


# ------------------------------------------------------------------------------
#
#
#                               VISUALIZATION UMAP
#
#
# ------------------------------------------------------------------------------

scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_1/umap_with_clusters.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="Reactive mature oligodendrocytes",
  "1"="Reactive myelinating oligodendrocytes",
  "2"="Activated vascular endothelial cells",
  "3"="Excitatory glutamatergic neurons",
  "4"="Disease-associated microglia (DAM)",
  "5"="Excitatory glutamatergic neurons",
  "6"="Infiltrating immune cells (lymphocyte/leukocyte)",
  "7"="Reactive astrocytes"
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 1
# -------------------------------------------------
# Convert cluster IDs to descriptive cell type labels
scanpy_data <- scanpy_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[1:length(cluster_labels)],
  cluster_labels
)

p_scanpy <- ggplot(scanpy_data,
                   aes(x = UMAP_1,
                       y = UMAP_2,
                       color = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  labs(
    title = "UMAP of Scanpy Cell Clusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_scanpy + NoLegend()

# -------------------------------------------------
#                      SPATIAL 1
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_1/umap_with_clusters.csv", show_col_types = FALSE)
cluster_vec <- scanpy_data$cluster

scanpy_data_spatial <- readr::read_csv("XeniumData_scanpy/scanpy_data_1/spatial_coords.csv", show_col_types = FALSE)
scanpy_data_spatial <- scanpy_data_spatial[,-4]
scanpy_data_spatial <- cbind(scanpy_data_spatial, cluster_vec)

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cluster_vec = as.character(cluster_vec),
    CellType = recode(cluster_vec, !!!cluster_labels),
    CellType = factor(
      CellType,
      levels = cluster_labels
    )
  )

p_scanpy_spatial <- ggplot(
    scanpy_data_spatial,
    aes(
      x = x,
      y = y,
      colour = CellType
    )) +
  
  geom_point( size = 0.45, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +      ## Xenium orientation
  labs( title = "Spatial Distribution of Scanpy Cell Types",
    x = "X Coordinate (µm)",
    y = "Y Coordinate (µm)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold" ),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right") +
  guides(colour = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1)))
p_scanpy_spatial
p_scanpy_spatial + NoLegend()


# -------------------------------------------------
#                      UMAP 2
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_2/umap_with_clusters.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="Disease-associated microglia (DAM)",
  "1"="Mature myelinating oligodendrocytes",
  "2"="Homeostatic microglia",
  "3"="Reactive astrocytes",
  "4"="Mature myelinating oligodendrocytes",
  "5"="Excitatory glutamatergic neurons",
  "6"="Vascular endothelial cells",
  "7"="FGFR3⁺ astrocytes",
  "8"="Oligodendrocyte precursor cells (OPCs / NG2 glia)",
  "9"="Excitatory glutamatergic neurons",
  "10"="LHX6⁺ inhibitory interneurons",
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# Convert cluster IDs to descriptive cell type labels
scanpy_data <- scanpy_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[1:length(cluster_labels)],
  cluster_labels
)

p_scanpy <- ggplot(scanpy_data,
                   aes(x = UMAP_1,
                       y = UMAP_2,
                       color = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  labs(
    title = "UMAP of Scanpy Cell Clusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_scanpy + NoLegend()

# -------------------------------------------------
#                      SPATIAL 2
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_2/umap_with_clusters.csv", show_col_types = FALSE)
cluster_vec <- scanpy_data$cluster

scanpy_data_spatial <- readr::read_csv("XeniumData_scanpy/scanpy_data_2/spatial_coords.csv", show_col_types = FALSE)
scanpy_data_spatial <- scanpy_data_spatial[,-4]
scanpy_data_spatial <- cbind(scanpy_data_spatial, cluster_vec)

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cluster_vec = as.character(cluster_vec),
    CellType = recode(cluster_vec, !!!cluster_labels),
    CellType = factor(
      CellType,
      levels = cluster_labels
    )
  )

p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    colour = CellType
  )) +
  
  geom_point( size = 0.45, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +      ## Xenium orientation
  labs( title = "Spatial Distribution of Scanpy Cell Types",
        x = "X Coordinate (µm)",
        y = "Y Coordinate (µm)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(
    face = "bold",
    hjust = 0.5,
    size = 16),
    axis.title = element_text(
      face = "bold" ),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right") +
  guides(colour = guide_legend(
    override.aes = list(
      size = 4,
      alpha = 1)))
p_scanpy_spatial
p_scanpy_spatial + NoLegend()


# -------------------------------------------------
#                      UMAP 3
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_3/umap_with_clusters.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="Reactive astrocytes",
  "1"="LHX6⁺ inhibitory interneurons (PV/SST lineage)",
  "2"="CCK⁺ excitatory glutamatergic neurons (HTR2A⁺ cortical neurons)",
  "3"="LAMP5⁺ upper-layer excitatory neurons",
  "4"="Excitatory glutamatergic projection neurons",
  "5"="Mature myelinating oligodendrocytes",
  "6"="Vascular endothelial cells",
  "7"="Oligodendrocyte precursor cells (OPCs / NG2 glia)",
  "8"="Activated microglia / Disease-associated microglia (DAM)",
  "9"="VIP⁺ / RELN⁺ inhibitory interneurons",
  "10"="LAMP5⁺ neurogliaform inhibitory interneurons"
)
cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# Convert cluster IDs to descriptive cell type labels
scanpy_data <- scanpy_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[1:length(cluster_labels)],
  cluster_labels
)

p_scanpy <- ggplot(scanpy_data,
                   aes(x = UMAP_1,
                       y = UMAP_2,
                       color = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  labs(
    title = "UMAP of Scanpy Cell Clusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_scanpy + NoLegend()

# -------------------------------------------------
#                      SPATIAL 3
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_3/umap_with_clusters.csv", show_col_types = FALSE)
cluster_vec <- scanpy_data$cluster

scanpy_data_spatial <- readr::read_csv("XeniumData_scanpy/scanpy_data_3/spatial_coords.csv", show_col_types = FALSE)
scanpy_data_spatial <- scanpy_data_spatial[,-4]
scanpy_data_spatial <- cbind(scanpy_data_spatial, cluster_vec)

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cluster_vec = as.character(cluster_vec),
    CellType = recode(cluster_vec, !!!cluster_labels),
    CellType = factor(
      CellType,
      levels = cluster_labels
    )
  )

p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    colour = CellType
  )) +
  
  geom_point( size = 0.45, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +      ## Xenium orientation
  labs( title = "Spatial Distribution of Scanpy Cell Types",
        x = "X Coordinate (µm)",
        y = "Y Coordinate (µm)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(
    face = "bold",
    hjust = 0.5,
    size = 16),
    axis.title = element_text(
      face = "bold" ),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right") +
  guides(colour = guide_legend(
    override.aes = list(
      size = 4,
      alpha = 1)))
p_scanpy_spatial
p_scanpy_spatial + NoLegend()

# -------------------------------------------------
#                      UMAP 3
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_4/umap_with_clusters.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="Vascular smooth muscle cells / Pericytes",
  "1"="Lonocyte-like epithelial cells (rare epithelial population; verify)",
  "2"="Activated / memory CD4⁺ T cells",
  "3"="Blood vascular endothelial cells",
  "4"="Activated macrophages (M2-like)",
  "5"="Plasma cells",
  "6"="Airway epithelial cells (secretory/luminal epithelial)",
  "7"="B lymphocytes",
  "8"="Alveolar type I (AT1) epithelial cells",
  "9"="Pulmonary capillary (aerocyte) endothelial cells",
  "10"="Alveolar type II (AT2) epithelial cells",
  "11"="Secretory / club epithelial cells",
  "12"="Mast cells"
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# Convert cluster IDs to descriptive cell type labels
scanpy_data <- scanpy_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[1:length(cluster_labels)],
  cluster_labels
)

p_scanpy <- ggplot(scanpy_data,
                   aes(x = UMAP_1,
                       y = UMAP_2,
                       color = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  labs(
    title = "UMAP of Scanpy Cell Clusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  )

p_scanpy + NoLegend()

# -------------------------------------------------
#                      SPATIAL 4
# -------------------------------------------------
scanpy_data <- readr::read_csv("XeniumData_scanpy/scanpy_data_4/umap_with_clusters.csv", show_col_types = FALSE)
cluster_vec <- scanpy_data$cluster

scanpy_data_spatial <- readr::read_csv("XeniumData_scanpy/scanpy_data_4/spatial_coords.csv", show_col_types = FALSE)
scanpy_data_spatial <- scanpy_data_spatial[,-4]
scanpy_data_spatial <- cbind(scanpy_data_spatial, cluster_vec)

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cluster_vec = as.character(cluster_vec),
    CellType = recode(cluster_vec, !!!cluster_labels),
    CellType = factor(
      CellType,
      levels = cluster_labels
    )
  )

p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    colour = CellType
  )) +
  
  geom_point( size = 0.45, alpha = 0.9) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +      ## Xenium orientation
  labs( title = "Spatial Distribution of Scanpy Cell Types",
        x = "X Coordinate (µm)",
        y = "Y Coordinate (µm)") +
  theme_classic(base_size = 15) +
  theme(plot.title = element_text(
    face = "bold",
    hjust = 0.5,
    size = 16),
    axis.title = element_text(
      face = "bold" ),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right") +
  guides(colour = guide_legend(
    override.aes = list(
      size = 4,
      alpha = 1)))
p_scanpy_spatial
p_scanpy_spatial + NoLegend()







