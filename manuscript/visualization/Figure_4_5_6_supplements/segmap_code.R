setwd("/Users/kxa643/Documents/xenium/segmap_R_codes")
getwd()

library(tidyverse)
library(FNN)
library(viridis)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(readr)
library(tibble)
library(tidyr)
library(purrr)
library(knitr)
library(Matrix)
library(spdep)
library(Seurat)

# ==========================================================
# Data 1
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_1/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="Mature oligodendrocytes",
  "1"="Excitatory neurons (Glutamatergic)",
  "2"="Astrocytes",
  "3"="Endothelial cells",
  "4"="PVALB interneurons",
  "5"="LAMP5/CUX2 excitatory neurons",
  "6"="OPCs (oligodendrocyte precursor cells)",
  "7"="Disease-associated microglia (DAM)"
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
segmap_data <- segmap_data %>%
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

p_segmap <- ggplot(segmap_data,
                   aes(x = umap1,
                       y = umap2,
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

p_segmap #+ NoLegend()

# -------------------------------------------------
#                      SPATIAL 1
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_1/data/spatial_clusters_final.csv", show_col_types = FALSE)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()

#---------------------
# Global Spatial Map
#---------------------

# -------------------------------------------------
#                      Local Diversity Analysis
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv(
  segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_1/data/spatial_clusters_final.csv", show_col_types = FALSE),
  show_col_types = FALSE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(cluster,
                     levels = unname(cluster_labels))
  )

# --------
# Compute Local Shannon Diversity
# --------

coords <- as.matrix(
  segmap_data_spatial[, c("x_centroid_px", "y_centroid_px")]
)

# Number of nearest neighbors
k_neighbors <- 50

knn_result <- get.knn(coords, k = k_neighbors)

# --------
# Shannon diversity function
# --------

shannon_diversity <- function(labels) {
  
  p <- table(labels) / length(labels)
  
  -sum(p * log(p))
}

# --------
# Calculate local diversity for each cell
# --------

local_diversity <- numeric(nrow(segmap_data_spatial))

for(i in seq_len(nrow(segmap_data_spatial))) {
  
  neighbor_ids <- knn_result$nn.index[i, ]
  
  neighbor_clusters <- segmap_data_spatial$cluster[neighbor_ids]
  
  local_diversity[i] <- shannon_diversity(neighbor_clusters)
}

segmap_data_spatial$local_diversity <- local_diversity

# --------
# Summary statistics
# --------

summary(segmap_data_spatial$local_diversity)

# --------
# Define diversity hotspots
# --------

threshold <- quantile(
  segmap_data_spatial$local_diversity,
  probs = 0.95,
  na.rm = TRUE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    hotspot = local_diversity >= threshold
  )

# --------
# Plot 1: Continuous Shannon Diversity Map
# --------
p_diversity <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px,
    color = local_diversity
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Shannon\nDiversity"
  ) +
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Cellular Diversity Map",
    subtitle = paste0(
      "k = ", k_neighbors,
      " nearest neighbors"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )
p_diversity

# --------
# Plot 2: Diversity Hotspots Only
# --------
p_hotspots <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px
  )
) +
  
  geom_point(
    color = "grey90",
    size = 0.25,
    alpha = 0.4
  ) +
  
  geom_point(
    data = subset(segmap_data_spatial, hotspot),
    color = "red",
    size = 0.6,
    alpha = 1
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Diversity Hotspots",
    subtitle = paste0(
      "Top 5% Shannon diversity regions"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
p_hotspots

# --------
# Plot 3: Distribution of Diversity Scores
# --------
p_hist <- ggplot(
  segmap_data_spatial,
  aes(local_diversity)
) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = threshold,
    linetype = 2,
    linewidth = 1
  ) +
  theme_classic(base_size = 15) +
  labs(
    title = "Distribution of Local Shannon Diversity",
    x = "Shannon Diversity",
    y = "Number of Cells"
  )
p_hist

# ==========================================================
# Data 2
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_2/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"  = "Myelinating_Oligodendrocytes",
  "1"  = "Activated_Microglia",
  "2"  = "Excitatory_Neurons_1 (IT/Corticocortical-like)",
  "3"  = "Reactive_Astrocytes",
  "4"  = "Endothelial_Cells",
  "5"  = "OPCs",
  "6"  = "Layer II/III Intratelencephalic-like",
  "7"  = "PV_Interneurons",
  "8"  = "HTR2A_Excitatory_Neurons",
  "9"  = "T_Cells",
  "10" = "Astrocyte_Progenitors"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 2
# -------------------------------------------------
# Convert cluster IDs to descriptive cell type labels
segmap_data <- segmap_data %>%
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

p_segmap <- ggplot(segmap_data,
                   aes(x = umap1,
                       y = umap2,
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

p_segmap# + NoLegend()

# -------------------------------------------------
#                      SPATIAL 2
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_2/data/spatial_clusters_final.csv", show_col_types = FALSE)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      Local Diversity Analysis
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv(
  segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_1/data/spatial_clusters_final.csv", show_col_types = FALSE),
  show_col_types = FALSE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(cluster,
                     levels = unname(cluster_labels))
  )

# --------
# Compute Local Shannon Diversity
# --------

coords <- as.matrix(
  segmap_data_spatial[, c("x_centroid_px", "y_centroid_px")]
)

# Number of nearest neighbors
k_neighbors <- 50

knn_result <- get.knn(coords, k = k_neighbors)

# --------
# Shannon diversity function
# --------

shannon_diversity <- function(labels) {
  
  p <- table(labels) / length(labels)
  
  -sum(p * log(p))
}

# --------
# Calculate local diversity for each cell
# --------

local_diversity <- numeric(nrow(segmap_data_spatial))

for(i in seq_len(nrow(segmap_data_spatial))) {
  
  neighbor_ids <- knn_result$nn.index[i, ]
  
  neighbor_clusters <- segmap_data_spatial$cluster[neighbor_ids]
  
  local_diversity[i] <- shannon_diversity(neighbor_clusters)
}

segmap_data_spatial$local_diversity <- local_diversity

# --------
# Summary statistics
# --------

summary(segmap_data_spatial$local_diversity)

# --------
# Define diversity hotspots
# --------

threshold <- quantile(
  segmap_data_spatial$local_diversity,
  probs = 0.95,
  na.rm = TRUE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    hotspot = local_diversity >= threshold
  )

# --------
# Plot 1: Continuous Shannon Diversity Map
# --------
p_diversity <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px,
    color = local_diversity
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Shannon\nDiversity"
  ) +
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Cellular Diversity Map",
    subtitle = paste0(
      "k = ", k_neighbors,
      " nearest neighbors"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )
p_diversity

# --------
# Plot 2: Diversity Hotspots Only
# --------
p_hotspots <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px
  )
) +
  
  geom_point(
    color = "grey90",
    size = 0.25,
    alpha = 0.4
  ) +
  
  geom_point(
    data = subset(segmap_data_spatial, hotspot),
    color = "red",
    size = 0.6,
    alpha = 1
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Diversity Hotspots",
    subtitle = paste0(
      "Top 5% Shannon diversity regions"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
p_hotspots

# --------
# Plot 3: Distribution of Diversity Scores
# --------
p_hist <- ggplot(
  segmap_data_spatial,
  aes(local_diversity)
) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = threshold,
    linetype = 2,
    linewidth = 1
  ) +
  theme_classic(base_size = 15) +
  labs(
    title = "Distribution of Local Shannon Diversity",
    x = "Shannon Diversity",
    y = "Number of Cells"
  )
p_hist



# ==========================================================
# Data 3
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_3/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"  = "Endothelial_Cells",
  "1"  = "Myelinating_Oligodendrocytes",
  "2"  = "OPCs",
  "3"  = "Activated_Microglia",
  "4"  = "RELN_Interneurons",
  "5"  = "VIP_Interneurons",
  "6"  = "RORB_Excitatory_Neurons",
  "7"  = "Astrocytes",
  "8"  = "CUX2_Excitatory_Neurons",
  "9"  = "SST_Interneurons",
  "10" = "PVALB_Interneurons"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 3
# -------------------------------------------------
# Convert cluster IDs to descriptive cell type labels
segmap_data <- segmap_data %>%
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

p_segmap <- ggplot(segmap_data,
                   aes(x = umap1,
                       y = umap2,
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

p_segmap #+ NoLegend()

# -------------------------------------------------
#                      SPATIAL 3
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_3/data/spatial_clusters_final.csv", show_col_types = FALSE)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      Local Diversity Analysis
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv(
  segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_1/data/spatial_clusters_final.csv", show_col_types = FALSE),
  show_col_types = FALSE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(cluster,
                     levels = unname(cluster_labels))
  )

# --------
# Compute Local Shannon Diversity
# --------

coords <- as.matrix(
  segmap_data_spatial[, c("x_centroid_px", "y_centroid_px")]
)

# Number of nearest neighbors
k_neighbors <- 50

knn_result <- get.knn(coords, k = k_neighbors)

# --------
# Shannon diversity function
# --------

shannon_diversity <- function(labels) {
  
  p <- table(labels) / length(labels)
  
  -sum(p * log(p))
}

# --------
# Calculate local diversity for each cell
# --------

local_diversity <- numeric(nrow(segmap_data_spatial))

for(i in seq_len(nrow(segmap_data_spatial))) {
  
  neighbor_ids <- knn_result$nn.index[i, ]
  
  neighbor_clusters <- segmap_data_spatial$cluster[neighbor_ids]
  
  local_diversity[i] <- shannon_diversity(neighbor_clusters)
}

segmap_data_spatial$local_diversity <- local_diversity

# --------
# Summary statistics
# --------

summary(segmap_data_spatial$local_diversity)

# --------
# Define diversity hotspots
# --------

threshold <- quantile(
  segmap_data_spatial$local_diversity,
  probs = 0.95,
  na.rm = TRUE
)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    hotspot = local_diversity >= threshold
  )

# --------
# Plot 1: Continuous Shannon Diversity Map
# --------
p_diversity <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px,
    color = local_diversity
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Shannon\nDiversity"
  ) +
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Cellular Diversity Map",
    subtitle = paste0(
      "k = ", k_neighbors,
      " nearest neighbors"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )
p_diversity

# --------
# Plot 2: Diversity Hotspots Only
# --------
p_hotspots <- ggplot(
  segmap_data_spatial,
  aes(
    x = x_centroid_px,
    y = y_centroid_px
  )
) +
  
  geom_point(
    color = "grey90",
    size = 0.25,
    alpha = 0.4
  ) +
  
  geom_point(
    data = subset(segmap_data_spatial, hotspot),
    color = "red",
    size = 0.6,
    alpha = 1
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  theme_classic(base_size = 15) +
  labs(
    title = "Local Diversity Hotspots",
    subtitle = paste0(
      "Top 5% Shannon diversity regions"
    ),
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )
p_hotspots

# --------
# Plot 3: Distribution of Diversity Scores
# --------
p_hist <- ggplot(
  segmap_data_spatial,
  aes(local_diversity)
) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "white"
  ) +
  geom_vline(
    xintercept = threshold,
    linetype = 2,
    linewidth = 1
  ) +
  theme_classic(base_size = 15) +
  labs(
    title = "Distribution of Local Shannon Diversity",
    x = "Shannon Diversity",
    y = "Number of Cells"
  )
p_hist


# ==========================================================
# Data 4
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_4/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"="B_cells",                     
  "1"="CD4_T_cells",                 
  "2"="Macrophages",                
  "3"="Epithelial tumor cells (EPCAM⁺)",          
  "4"="Epithelial tumor cells (EGFR⁺ secretory subtype)",         
  "5"="Plasma_cells",                
  "6"="Mast_cells",                  
  "7"="Aerocyte_endothelial",        
  "8"="Endothelial_cells",           
  "9"="Secretory_epithelial",        
  "10"="Pericyte_SMC",                
  "11"="Monocytes",                   
  "12"="Goblet_mucous_cells"          
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 4
# -------------------------------------------------

# Convert cluster IDs to descriptive cell type labels
segmap_data <- segmap_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  cluster_labels
)

# Create Seurat object
dummy_counts <- Matrix::Matrix(
  0,
  nrow = 1,
  ncol = nrow(segmap_data),
  sparse = TRUE
)

colnames(dummy_counts) <- segmap_data$cell_id
rownames(dummy_counts) <- "dummy"

segmap <- CreateSeuratObject(dummy_counts)

# Add metadata
segmap$cluster <- segmap_data$cluster

# Add UMAP coordinates
umap_embeddings <- as.matrix(
  segmap_data[, c("umap1", "umap2")]
)

rownames(umap_embeddings) <- segmap_data$cell_id
colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")

segmap[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings,
  key = "UMAP_",
  assay = DefaultAssay(segmap)
)


# Set identities
Idents(segmap) <- "cluster"

# Plot
p_segmap_spatial <- DimPlot(
  segmap,
  reduction = "umap",
  cols = celltype_colors,
  pt.size = 0.35,
  raster=FALSE
) +
  ggtitle("UMAP of SegMap Cell Clusters") +
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
  )

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      SPATIAL 4
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_4/data/spatial_clusters_final.csv", show_col_types = FALSE)
segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()


# ==========================================================
# SEGMAP 5
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_5/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"  = "B_Cells",
  "1"  = "CD8_T_Cells",
  "2"  = "Monocytes",
  "3"  = "CD4_T_Cells",
  "4"  = "iCAFs",
  "5"  = "Plasma_Cells",
  "6"  = "Fibroblasts",
  "7"  = "Macrophages",
  "8"  = "Smooth_Muscle_Pericytes",
  "9"  = "Tumor_Epithelial_Cells",
  "10" = "Venous_Endothelial_Cells",
  "11" = "Mast_Cells",
  "12" = "AT1_Cells",
  "13" = "AT2_Cells: Secretory Epithelial Cells",
  "14" = "Pericytes",
  "15" = "Alveolar_Macrophages",
  "16" = "Cycling_Tumor_Cells",
  "17" = "Arterial_Endothelial_Cells",
  "18" = "Capillary_Endothelial_Cells",
  "19" = "AT2_like_Tumor_Cells",
  "20" = "Secretory_Tumor_Cells"
)


cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 5
# -------------------------------------------------

# Convert cluster IDs to descriptive cell type labels
segmap_data <- segmap_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  cluster_labels
)

# Create Seurat object
dummy_counts <- Matrix::Matrix(
  0,
  nrow = 1,
  ncol = nrow(segmap_data),
  sparse = TRUE
)

colnames(dummy_counts) <- segmap_data$cell_id
rownames(dummy_counts) <- "dummy"

segmap <- CreateSeuratObject(dummy_counts)

# Add metadata
segmap$cluster <- segmap_data$cluster

# Add UMAP coordinates
umap_embeddings <- as.matrix(
  segmap_data[, c("umap1", "umap2")]
)

rownames(umap_embeddings) <- segmap_data$cell_id
colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")

segmap[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings,
  key = "UMAP_",
  assay = DefaultAssay(segmap)
)


# Set identities
Idents(segmap) <- "cluster"

# Plot
p_segmap_spatial <- DimPlot(
  segmap,
  reduction = "umap",
  cols = celltype_colors,
  pt.size = 0.35,
  raster=FALSE
) +
  ggtitle("UMAP of SegMap Cell Clusters") +
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
  )

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      SPATIAL 5
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_5/data/spatial_clusters_final.csv", show_col_types = FALSE)


segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      ZOOM
# -------------------------------------------------

# -------------
# CONVERT CLUSTERS TO CELL-TYPE LABELS
# -------------

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

# -------------
# NAMED CELL-TYPE COLORS
# -------------

celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

# -------------
# DEFINE ROI
# -------------
roi_xmin <- 18000
roi_xmax <- 28000

roi_ymin <- 12000
roi_ymax <- 22000

# -------------
# EXTRACT ROI
# -------------

roi_df <- segmap_data_spatial %>%
  filter(
    x_centroid_px >= roi_xmin,
    x_centroid_px <= roi_xmax,
    y_centroid_px >= roi_ymin,
    y_centroid_px <= roi_ymax
  )

cat("Cells in ROI:", nrow(roi_df), "\n")

# -------------
# PANEL A: WHOLE TISSUE MAP
# -------------

p_full <-
  
  ggplot(
    segmap_data_spatial,
    aes(
      x = x_centroid_px,
      y = y_centroid_px,
      colour = cluster
    )
  ) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85
  ) +
  
  annotate(
    "rect",
    xmin = roi_xmin,
    xmax = roi_xmax,
    ymin = roi_ymin,
    ymax = roi_ymax,
    colour = "red",
    fill = NA,
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  
  labs(
    title = "Whole Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "none"
  )
p_full

# -------------
# PANEL B: MAGNIFIED ROI
# -------------

p_zoom <-
  
  ggplot(
    roi_df,
    aes(
      x = x_centroid_px,
      y = y_centroid_px,
      colour = cluster
    )
  ) +
  
  geom_point(
    size = 1.3,
    alpha = 0.95
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  
  labs(
    title = "Magnified ROI",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key.height = unit(
      0.5,
      "cm"
    ),
    legend.position = "none"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )
p_zoom

# -------------
# COMBINE PANELS
# -------------

final_plot <-
  
  p_full +
  p_zoom +
  plot_layout(
    widths = c(1.3, 1)
  )
final_plot


# ==========================================================
# SEGMAP 6
# ==========================================================
segmap_data <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_6/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_labels <- c(
  "0"  = "Mito_High_Cells",
  "1"  = "Reactive_Oligodendrocytes",
  "2"  = "Perivascular_Fibroblasts",
  "3"  = "Astrocytes",
  "4"  = "Excitatory_Neurons",
  "5"  = "Excitatory_Neurons_MEF2C",
  "6"  = "Myelinating_Oligodendrocytes",
  "7"  = "Perivascular_Macrophages",
  "8"  = "RORB_Excitatory_Neurons",
  "9"  = "Endothelial_Cells",
  "10" = "LAMP5_CUX2_Neurons",
  "11" = "OPCs",
  "12" = "Activated_Microglia",
  "13" = "Reactive_Glia",
  "14" = "Neurogliaform_Related_Neurons",
  "15" = "PVALB_Interneurons"
)

cluster_colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#FF0000", "#FFE4E1",
  "#8B0000", "#E377C2", "#00FFFF", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#00FF7F", "#551A8B", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#636363", "#E6550D",
  "#0000FF", "#CD7054", "#66A61E", "#E6AB02", "#A6761D",
  "#1B9E77", "#D95F02", "#FF00FF", "#000000")

# -------------------------------------------------
#                      UMAP 6
# -------------------------------------------------

# Convert cluster IDs to descriptive cell type labels
segmap_data <- segmap_data %>%
  mutate(
    cluster = factor(
      as.character(cluster),
      levels = names(cluster_labels),
      labels = cluster_labels
    )
  )

# Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  cluster_labels
)

# Create Seurat object
dummy_counts <- Matrix::Matrix(
  0,
  nrow = 1,
  ncol = nrow(segmap_data),
  sparse = TRUE
)

colnames(dummy_counts) <- segmap_data$cell_id
rownames(dummy_counts) <- "dummy"

segmap <- CreateSeuratObject(dummy_counts)

# Add metadata
segmap$cluster <- segmap_data$cluster

# Add UMAP coordinates
umap_embeddings <- as.matrix(
  segmap_data[, c("umap1", "umap2")]
)

rownames(umap_embeddings) <- segmap_data$cell_id
colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")

segmap[["umap"]] <- CreateDimReducObject(
  embeddings = umap_embeddings,
  key = "UMAP_",
  assay = DefaultAssay(segmap)
)


# Set identities
Idents(segmap) <- "cluster"

# Plot
p_segmap_spatial <- DimPlot(
  segmap,
  reduction = "umap",
  cols = celltype_colors,
  pt.size = 0.35,
  raster=FALSE
) +
  ggtitle("UMAP of SegMap Cell Clusters") +
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
  )

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      SPATIAL 6
# -------------------------------------------------
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/cellpose_pipeline_outputs_6/data/spatial_clusters_final.csv", show_col_types = FALSE)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

## Named color vector
celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85) +
  
  scale_color_manual(values = celltype_colors,
                     drop = FALSE,
                     name = "Cell type") +
  
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "SegMap Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "right"
  ) +
  
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))

p_segmap_spatial
p_segmap_spatial + NoLegend()

# -------------------------------------------------
#                      ZOOM
# -------------------------------------------------

# -------------
# CONVERT CLUSTERS TO CELL-TYPE LABELS
# -------------

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_labels),
    cluster = factor(
      cluster,
      levels = unname(cluster_labels)
    )
  )

# -------------
# NAMED CELL-TYPE COLORS
# -------------

celltype_colors <- setNames(
  cluster_colors[seq_along(cluster_labels)],
  unname(cluster_labels)
)

# -------------
# DEFINE ROI
# -------------
roi_xmin <- 18000
roi_xmax <- 28000

roi_ymin <- 12000
roi_ymax <- 22000

# -------------
# EXTRACT ROI
# -------------

roi_df <- segmap_data_spatial %>%
  filter(
    x_centroid_px >= roi_xmin,
    x_centroid_px <= roi_xmax,
    y_centroid_px >= roi_ymin,
    y_centroid_px <= roi_ymax
  )

cat("Cells in ROI:", nrow(roi_df), "\n")

# -------------
# PANEL A: WHOLE TISSUE MAP
# -------------

p_full <-
  
  ggplot(
    segmap_data_spatial,
    aes(
      x = x_centroid_px,
      y = y_centroid_px,
      colour = cluster
    )
  ) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85
  ) +
  
  annotate(
    "rect",
    xmin = roi_xmin,
    xmax = roi_xmax,
    ymin = roi_ymin,
    ymax = roi_ymax,
    colour = "red",
    fill = NA,
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  
  labs(
    title = "Whole Spatial Cell-Type Map",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16),
    axis.title = element_text(
      face = "bold"),
    legend.title = element_text(
      face = "bold"),
    legend.text = element_text(
      size = 10),
    legend.key.height = unit(0.5, "cm"),
    legend.position = "none"
  ) +
  guides(colour = guide_legend(
    override.aes = list(
      size = 3,
      alpha = 1)))
p_full

# -------------
# PANEL B: MAGNIFIED ROI
# -------------

p_zoom <-
  
  ggplot(
    roi_df,
    aes(
      x = x_centroid_px,
      y = y_centroid_px,
      colour = cluster
    )
  ) +
  
  geom_point(
    size = 1.3,
    alpha = 0.95
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell type"
  ) +
  
  coord_equal() +
  scale_y_reverse() +
  
  labs(
    title = "Magnified ROI",
    x = "X coordinate (pixels)",
    y = "Y coordinate (pixels)"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key.height = unit(
      0.5,
      "cm"
    ),
    legend.position = "none"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )
p_zoom

# -------------
# COMBINE PANELS
# -------------

final_plot <-
  
  p_full +
  p_zoom +
  plot_layout(
    widths = c(1.3, 1)
  )
final_plot




