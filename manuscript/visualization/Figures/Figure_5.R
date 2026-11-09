
setwd("/Users/kxa643/Documents/CWRU/xenium")
getwd()

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid)


# ==============================================================================
# SCANPY
# ==============================================================================
load("CLUSTERS_MAXIMUM/Scanpy_High_Cluster.rdata")
# --------------------------------
# Harmonized cluster annotation
# --------------------------------

cluster_annotation_scanpy_harmonized <- c(
  "0" = "Oligodendrocyte_Mature",
  "1" = "Excitatory_IT_HS3ST2",
  "2" = "Interneuron_PVALB",
  "3" = "Astrocyte_Homeostatic",
  "4" = "Excitatory_LAMP5_CUX2",
  "5" = "Microglia_DAM_SPP1",
  "6" = "Endothelial_FLT1",
  "7" = "OPC_PDGFRA",
  "8" = "Interneuron_LAMP5_KIT"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Oligodendrocyte_Mature" = "#1F77B4",
  "Excitatory_IT_HS3ST2"   = "#FF7F0E",
  "Interneuron_PVALB"      = "#2CA02C",
  "Astrocyte_Homeostatic"  = "#D62728",
  "Excitatory_LAMP5_CUX2"  = "#9467BD",
  "Microglia_DAM_SPP1"     = "#8C564B",
  "Endothelial_FLT1"       = "#E377C2",
  "OPC_PDGFRA"             = "#17BECF",
  "Interneuron_LAMP5_KIT"  = "#E6550D"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

scanpy_data <- as.data.frame(
  Embeddings(
    Scanpy_High_Cluster,
    reduction = "umap"
  )
)

colnames(scanpy_data)[1:2] <- c("UMAP_1", "UMAP_2")

scanpy_data$cluster <- as.character(
  Idents(Scanpy_High_Cluster)
)

# --------------------------------
# Map cluster -> cell type
# --------------------------------

scanpy_data$cell_type <- cluster_annotation_scanpy_harmonized[
  scanpy_data$cluster
]

scanpy_data$cell_type <- factor(
  scanpy_data$cell_type,
  levels = names(celltype_colors)
)

# --------------------------------
# Compute cluster centers
# --------------------------------

cluster_centers <- scanpy_data %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# --------------------------------
# UMAP plot
# --------------------------------

p_scanpy <- ggplot(
  scanpy_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = cell_type
  )
) +
  
  geom_point(
    size = 0.25,
    alpha = 0.8
  ) +
  
  #  geom_label_repel(
  #    data = cluster_centers,
  #    aes(
  #      x = UMAP_1,
  #      y = UMAP_2,
  #      label = cell_type
  #    ),
  #    inherit.aes = FALSE,
  #    fontface = "bold",
  #    size = 4,
  #    color = "black",
  #    fill = "white",
  #    box.padding = 0.5,
  #    point.padding = 0.3,
  #    segment.color = "grey50",
  #    max.overlaps = Inf,
  #    show.legend = FALSE
  #  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  
  labs(
    title = "Scanpy UMAP (Harmonized Cell Types)",
    x = "umap 1",
    y = "umap 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.key.height = unit(
      0.5,
      "cm"
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

# --------------------------------
# Display plot
# --------------------------------

p_scanpy

# --------------------------------
# Optional: hide legend
# --------------------------------

p_scanpy +
  theme(
    legend.position = "none"
  )

#===================================================
# spatial visualization
#===================================================

# -----------------------------
# Prepare metadata
# -----------------------------
spatial_data <- Scanpy_High_Cluster@meta.data

spatial_data$cell_type <- factor(
  cluster_annotation_scanpy_harmonized[
    as.character(spatial_data$leiden)
  ],
  levels = names(celltype_colors)
)

# -----------------------------
# Cell-type centroids for labels
# -----------------------------

cluster_centers <- spatial_data %>%
  group_by(cell_type) %>%
  summarize(
    x_centroid = median(x_centroid, na.rm = TRUE),
    y_centroid = median(y_centroid, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Spatial plot
# -----------------------------

p_spatial_scanpy <- ggplot(
  spatial_data,
  aes(
    x = x_centroid,
    y = y_centroid,
    color = cell_type
  )
) +
  
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  
  #  geom_label_repel(
  #    data = cluster_centers,
  #    aes(
  #      x = x_centroid,
  #      y = y_centroid,
  #      label = cell_type
  #    ),
  #    inherit.aes = FALSE,
  #    size = 4,
  #    fontface = "bold",
  #    color = "black",
  #    fill = "white",
  #    box.padding = 0.6,
  #    point.padding = 0.4,
  #    segment.color = "grey50",
  #    segment.size = 0.4,
  #    show.legend = FALSE,
  #    max.overlaps = Inf
  #  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE
  ) +
  
  scale_y_reverse() +
  
  coord_fixed() +
  
  labs(
    title = "Spatial Distribution of Harmonized Cell Types",
    x = NULL,
    y = NULL,
    color = "Cell Type"
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    axis.line = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

p_spatial_scanpy#+
#theme(
#  legend.position = "none"
#)


# ==============================================================================
# SEURAT
# ==============================================================================
load("CLUSTERS_MAXIMUM/Seurat_High_Cluster.rdata")

# --------------------------------
# Harmonized cluster annotation
# --------------------------------

cluster_annotation_seurat_harmonized <- c(
  "1"  = "Oligodendrocyte_Mature",
  "2"  = "Oligodendrocyte_Reactive_SPP1",
  "3"  = "Endothelial_FLT1",
  "4"  = "Excitatory_IT_HS3ST2",
  "5"  = "Astrocyte_Homeostatic",
  "6"  = "Microglia_DAM_SPP1",
  "7"  = "Interneuron_PVALB",
  "8"  = "OPC_PDGFRA",
  "9"  = "Oligodendrocyte_Reactive_CRYAB",
  "10" = "Excitatory_RORB",
  "11" = "Excitatory_Cycling",
  "12" = "Excitatory_LAMP5_CUX2",
  "13" = "Astrocyte_Reactive_GFAP",
  "14" = "Interneuron_SST"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Oligodendrocyte_Mature"         = "#1F77B4",
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Endothelial_FLT1"              = "#E377C2",
  "Excitatory_IT_HS3ST2"          = "#FF7F0E",
  "Astrocyte_Homeostatic"         = "#D62728",
  "Microglia_DAM_SPP1"            = "#8C564B",
  "Interneuron_PVALB"             = "#2CA02C",
  "OPC_PDGFRA"                    = "#17BECF",
  "Oligodendrocyte_Reactive_CRYAB" = "#FFE4E1",
  "Excitatory_RORB"               = "#551A8B",
  "Excitatory_Cycling"            = "#00FFFF",
  "Excitatory_LAMP5_CUX2"         = "#9467BD",
  "Astrocyte_Reactive_GFAP"       = "#00FF7F",
  "Interneuron_SST"               = "#8B0000"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

seurat_data <- as.data.frame(
  Embeddings(
    Seurat_High_Cluster,
    reduction = "umap"
  )
)

colnames(seurat_data)[1:2] <- c("UMAP_1", "UMAP_2")

seurat_data$cluster <- as.character(
  Idents(Seurat_High_Cluster)
)

# --------------------------------
# Map cluster -> cell type
# --------------------------------

seurat_data$cell_type <- cluster_annotation_seurat_harmonized[
  seurat_data$cluster
]

seurat_data$cell_type <- factor(
  seurat_data$cell_type,
  levels = names(celltype_colors)
)

# --------------------------------
# Compute cluster centers
# --------------------------------

cluster_centers <- seurat_data %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# --------------------------------
# UMAP plot
# --------------------------------

p_seurat <- ggplot(
  seurat_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = cell_type
  )
) +
  
  geom_point(
    size = 0.25,
    alpha = 0.8
  ) +
  
  #  geom_label_repel(
  #    data = cluster_centers,
  #    aes(
  #      x = UMAP_1,
  #      y = UMAP_2,
  #      label = cell_type
  #    ),
  #    inherit.aes = FALSE,
  #    fontface = "bold",
  #    size = 4,
  #    color = "black",
  #    fill = "white",
  #    box.padding = 0.5,
  #    point.padding = 0.3,
  #    segment.color = "grey50",
  #    max.overlaps = Inf,
  #    show.legend = FALSE
  #  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  
  labs(
    title = "Seurat UMAP (Harmonized Cell Types)",
    x = "umap 1",
    y = "umap 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.key.height = unit(
      0.5,
      "cm"
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

# --------------------------------
# Display plot
# --------------------------------

p_seurat

# --------------------------------
# Optional: hide legend
# --------------------------------

p_seurat +
  theme(
    legend.position = "none"
  )

#===================================================
# spatial visualization
#===================================================
# Get coordinates and cluster IDs
coords <- GetTissueCoordinates(Seurat_High_Cluster)
coords$cluster <- Idents(Seurat_High_Cluster)

# Convert cluster to factor with desired order
coords$cluster <- factor(
  coords$cluster,
  levels = names(cluster_annotation_seurat_harmonized),
  labels = cluster_annotation_seurat_harmonized
)

p_seurat_spatial <- ggplot(coords, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_color_manual(
    values = setNames(celltype_colors, cluster_annotation_seurat_harmonized),
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


# ==============================================================================
# BAYSOR
# ==============================================================================
load("CLUSTERS_MAXIMUM/Baysor_High_Cluster.rdata")

cluster_annotation_baysor_harmonized <- c(
  "0"  = "Oligodendrocyte_Mature",
  "1"  = "Oligodendrocyte_Reactive_SPP1",
  "2"  = "Excitatory_RORB",
  "3"  = "Astrocyte_Homeostatic",
  "4"  = "Endothelial_FLT1",
  "5"  = "Interneuron_PVALB",
  "6"  = "Astrocyte_Reactive_GFAP",
  "7"  = "Microglia_DAM_SPP1",
  "8"  = "Excitatory_LAMP5_CUX2",
  "9"  = "OPC_PDGFRA",
  "10" = "Interneuron_VIP_PROX1",
  "11" = "Oligodendrocyte_MOBP_SULF1",
  "12" = "SPP1_Unclassified"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Oligodendrocyte_Mature"        = "#1F77B4",
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Excitatory_RORB"               = "#551A8B",
  "Astrocyte_Homeostatic"         = "#D62728",
  "Endothelial_FLT1"              = "#E377C2",
  "Interneuron_PVALB"             = "#2CA02C",
  "Astrocyte_Reactive_GFAP"       = "#00FF7F",
  "Microglia_DAM_SPP1"            = "#8C564B",
  "Excitatory_LAMP5_CUX2"         = "#9467BD",
  "OPC_PDGFRA"                    = "#17BECF",
  "Interneuron_VIP_PROX1"         = "#BCBD22",
  "Oligodendrocyte_MOBP_SULF1"    = "#CD7054",
  "SPP1_Unclassified"             = "#000000"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

baysor_data <- as.data.frame(
  Embeddings(
    Baysor_High_Cluster,
    reduction = "umap"
  )
)

colnames(baysor_data)[1:2] <- c("UMAP_1", "UMAP_2")

baysor_data$cluster <- as.character(
  Idents(Baysor_High_Cluster)
)

# --------------------------------
# Map cluster -> cell type
# --------------------------------

baysor_data$cell_type <- cluster_annotation_baysor_harmonized[
  baysor_data$cluster
]

baysor_data$cell_type <- factor(
  baysor_data$cell_type,
  levels = names(celltype_colors)
)

# --------------------------------
# Compute cluster centers
# --------------------------------

cluster_centers <- baysor_data %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# --------------------------------
# UMAP plot
# --------------------------------

p_baysor <- ggplot(
  baysor_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = cell_type
  )
) +
  
  geom_point(
    size = 0.25,
    alpha = 0.8
  ) +
  
  #  geom_label_repel(
  #    data = cluster_centers,
  #    aes(
  #      x = UMAP_1,
  #      y = UMAP_2,
  #      label = cell_type
  #    ),
  #    inherit.aes = FALSE,
  #    fontface = "bold",
  #    size = 4,
  #    color = "black",
  #    fill = "white",
  #    box.padding = 0.5,
  #    point.padding = 0.3,
  #    segment.color = "grey50",
  #    max.overlaps = Inf,
  #    show.legend = FALSE
  #  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  
  labs(
    title = "Baysor UMAP (Harmonized Cell Types)",
    x = "umap 1",
    y = "umap 2"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 18
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    legend.position = "right",
    
    legend.key.height = unit(
      0.5,
      "cm"
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

# --------------------------------
# Display plot
# --------------------------------

p_baysor

# --------------------------------
# Optional: hide legend
# --------------------------------

p_baysor +
  theme(
    legend.position = "none"
  )

#===================================================
# spatial visualization
#===================================================

# -------------------------------------------------
#                      SPATIAL 1
# -------------------------------------------------
celltype_colors <- setNames(
  celltype_colors,
  unname(cluster_annotation_baysor_harmonized)
)

## Metadata
spatial.df <- Baysor_High_Cluster@meta.data
spatial.df$cell <- rownames(spatial.df)
spatial.df$cluster <- as.character(spatial.df$seurat_clusters)

## Add cell-type labels
spatial.df$cell_type <- cluster_annotation_baysor_harmonized[spatial.df$cluster]

## Remove cells without annotations (optional)
spatial.df <- spatial.df[!is.na(spatial.df$cell_type), ]

## Keep legend in desired order
spatial.df$cell_type <- factor(
  spatial.df$cell_type,
  levels = unname(cluster_annotation_baysor_harmonized)
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
  geom_point(size = 0.25, alpha = 0.9) +
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

# ==============================================================================
# SEGMAP
# ==============================================================================

segmap_data <- readr::read_csv("CLUSTERS_MAXIMUM/cellpose_pipeline_outputs_1_iter15/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_annotation_segmap_harmonized <- c(
  "0"  = "Astrocyte_Reactive_GFAP",
  "1"  = "Astrocyte_Homeostatic",
  "2"  = "Oligodendrocyte_Reactive_SPP1",
  "3"  = "Oligodendrocyte_MOBP_Stressed",
  "4"  = "Excitatory_LAMP5_CUX2",
  "5"  = "Interneuron_VIP_PROX1",
  "6"  = "Interneuron_LAMP5_KIT",
  "7"  = "Oligodendrocyte_Mature",
  "8"  = "Excitatory_IT_HS3ST2",
  "9"  = "Excitatory_IT_MEF2C",
  "10" = "Excitatory_NR4A2_TRHDE",
  "11" = "OPC_PDGFRA",
  "12" = "Pericyte_ABCC9",
  "13" = "Excitatory_RORB",
  "14" = "Interneuron_PVALB",
  "15" = "Microglia_DAM_SPP1",
  "16" = "Endothelial_FLT1"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Astrocyte_Reactive_GFAP"       = "#00FF7F",
  "Astrocyte_Homeostatic"         = "#D62728",
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Oligodendrocyte_MOBP_Stressed" = "#FF00FF",
  "Excitatory_LAMP5_CUX2"         = "#9467BD",
  "Interneuron_VIP_PROX1"         = "#BCBD22",
  "Interneuron_LAMP5_KIT"         = "#E6550D",
  "Oligodendrocyte_Mature"        = "#1F77B4",
  "Excitatory_IT_HS3ST2"          = "#FF7F0E",
  "Excitatory_IT_MEF2C"           = "#31A354",
  "Excitatory_NR4A2_TRHDE"        = "#E6AB02",
  "OPC_PDGFRA"                    = "#17BECF",
  "Pericyte_ABCC9"                = "#7B4173",
  "Excitatory_RORB"               = "#551A8B",
  "Interneuron_PVALB"             = "#2CA02C",
  "Microglia_DAM_SPP1"            = "#8C564B",
  "Endothelial_FLT1"              = "#E377C2"
)

#===================================================
# Convert cluster IDs to harmonized cell types
#===================================================

segmap_data <- segmap_data %>%
  mutate(
    cluster = as.character(cluster),
    
    cell_type = factor(
      cluster_annotation_segmap_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

#===================================================
# UMAP Plot
#===================================================

p_segmap <- ggplot(
  segmap_data,
  aes(
    x = umap1,
    y = umap2,
    color = cell_type
  )
) +
  
  geom_point(
    size = 0.25,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  
  labs(
    title = "UMAP of SegMap Cell Types",
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
    
    legend.position = "right"
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size =4,
        alpha = 1
      )
    )
  )

p_segmap #+ NoLegend()

#===================================================
# spatial visualization
#===================================================
segmap_data_spatial <- readr::read_csv("CLUSTERS_MAXIMUM/cellpose_pipeline_outputs_1_iter15/data/spatial_clusters_final.csv", show_col_types = FALSE)

segmap_data_spatial <- segmap_data_spatial %>%
  mutate(
    cluster = as.character(cluster),
    cluster = recode(cluster, !!!cluster_annotation_segmap_harmonized),
    cluster = factor(
      cluster,
      levels = unname(cluster_annotation_segmap_harmonized)
    )
  )

## Named color vector
celltype_colors <- setNames(
  celltype_colors[seq_along(cluster_annotation_segmap_harmonized)],
  unname(cluster_annotation_segmap_harmonized)
)

## Spatial plot
p_segmap_spatial <-
  
  ggplot(segmap_data_spatial,
         aes( x = x_centroid_px,
              y = y_centroid_px,
              colour = cluster)) +
  
  geom_point(
    size = 0.25,
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