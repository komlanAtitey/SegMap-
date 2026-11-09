library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)

# ------------------------------------------------------------------------------
#
#
# DATA 1
#
#
# ------------------------------------------------------------------------------


# ==============================================================================
# SEGMAP
# ==============================================================================

segmap_data <- readr::read_csv("XeniumData_segmap/segmapOutput_1/cellpose_pipeline_outputs_1_iter15/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_annotation_segmap_harmonized <- c(
  "0" = "Astrocyte_Reactive_GFAP",
  "1" = "Oligodendrocyte_Reactive_SPP1",
  "2" = "Excitatory_IT_MEF2C",
  "3" = "Endothelial_FLT1",
  "4" = "OPC_PDGFRA",
  "5" = "Interneuron_PVALB",
  "6" = "Excitatory_LAMP5_CUX2",
  "7" = "Microglia_DAM_SPP1"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Astrocyte_Reactive_GFAP"       = "#00FF7F",
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Excitatory_IT_MEF2C"           = "#31A354",
  "Endothelial_FLT1"              = "#E377C2",
  "OPC_PDGFRA"                    = "#17BECF",
  "Interneuron_PVALB"             = "#2CA02C",
  "Excitatory_LAMP5_CUX2"         = "#9467BD",
  "Microglia_DAM_SPP1"            = "#8C564B"
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
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/segmapOutput_1/cellpose_pipeline_outputs_1_iter15/data/spatial_clusters_final.csv", show_col_types = FALSE)
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

# ==============================================================================
#
#
# DATA 2
#
#
# ==============================================================================

segmap_data <- readr::read_csv("XeniumData_segmap/segmapOutput_2/cellpose_pipeline_outputs_1_iter15/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_annotation_segmap_harmonized <- c(
  "0"  = "Astrocyte_Reactive_APOE",
  "1"  = "Microglia_DAM_ITGAX",
  "2"  = "OPC_PDGFRA",
  "3"  = "Oligodendrocyte_Mature",
  "4"  = "Endothelial_FLT1",
  "5"  = "Excitatory_C1QL3_EPHA4",
  "6"  = "Excitatory_RORB_HS3ST2",
  "7"  = "VIP_PROX1_Interneuron",
  "8"  = "SST_LHX6_Interneuron",
  "9"  = "PVALB_TAC1_Interneuron",
  "10" = "Excitatory_TSHZ2_HTR2C"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Astrocyte_Reactive_APOE"     = "#3182BD",
  "Microglia_DAM_ITGAX"         = "#CD7054",
  "OPC_PDGFRA"                  = "#31A354",
  "Oligodendrocyte_Mature"      = "#2CA02C",
  "Endothelial_FLT1"            = "#0000FF",
  "Excitatory_C1QL3_EPHA4"      = "#637939",
  "Excitatory_RORB_HS3ST2"      = "#551A8B",
  "VIP_PROX1_Interneuron"       = "#00FF7F",
  "SST_LHX6_Interneuron"        = "#00FFFF",
  "PVALB_TAC1_Interneuron"      = "#66A61E",
  "Excitatory_TSHZ2_HTR2C"      = "#17BECF"
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
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/segmapOutput_2/cellpose_pipeline_outputs_1_iter15/data/spatial_clusters_final.csv", show_col_types = FALSE)
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

# ==============================================================================
#
#
# DATA 3
#
#
# ==============================================================================

segmap_data <- readr::read_csv("XeniumData_segmap/segmapOutput_3/cellpose_pipeline_outputs_1_iter15/data/umap_clusters_final.csv", show_col_types = FALSE)


cluster_annotation_segmap_harmonized <- c(
  "0"  = "Astrocyte_Reactive",
  "1"  = "Microglia_DAM",
  "2"  = "OPC_PDGFRA",
  "3"  = "Oligodendrocyte_Mature",
  "4"  = "Endothelial_FLT1",
  "5"  = "Excitatory_LAMP5_C1QL3",
  "6"  = "Excitatory_RORB_HS3ST2",
  "7"  = "Interneuron_VIP_PROX1",
  "8"  = "Interneuron_SST_LHX6",
  "9"  = "Interneuron_PVALB_LHX6",
  "10" = "Excitatory_TSHZ2_CRYM"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Astrocyte_Reactive"       = "#3182BD",
  "Microglia_DAM"            = "#FFE4E1",
  "OPC_PDGFRA"               = "#31A354",
  "Oligodendrocyte_Mature"   = "#2CA02C",
  "Endothelial_FLT1"         = "#0000FF",
  "Excitatory_LAMP5_C1QL3"   = "#E6550D",
  "Excitatory_RORB_HS3ST2"   = "#551A8B",
  "Interneuron_VIP_PROX1"    = "#8B0000",
  "Interneuron_SST_LHX6"     = "#E377C2",
  "Interneuron_PVALB_LHX6"   = "#756BB1",
  "Excitatory_TSHZ2_CRYM"    = "#1B9E77"
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
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/segmapOutput_3/cellpose_pipeline_outputs_1_iter15/data/spatial_clusters_final.csv", show_col_types = FALSE)
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

# ==============================================================================
#
#
# DATA 4
#
#
# ==============================================================================

segmap_data <- readr::read_csv("XeniumData_segmap/segmapOutput_4/cellpose_pipeline_outputs_1_iter15/data/umap_clusters_final.csv", show_col_types = FALSE)

cluster_annotation_segmap_harmonized <- c(
  "0"  = "CD4_Tcell_Activated",
  "1"  = "Aerocyte_Endothelial",
  "2"  = "Bcell_Naive",
  "3"  = "Plasma_Cell_MZB1",
  "4"  = "Fibroblast_Matrix",
  "5"  = "Vascular_Endothelial",
  "6"  = "Pericyte_SMC",
  "7"  = "Mast_Cell",
  "8"  = "Alveolar_Macrophage",
  "9"  = "Club_Ciliated_Transitional",
  "10" = "Epithelial_Progenitor",
  "11" = "Goblet_MUC5B",
  "12" = "AT2_Transitional"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "CD4_Tcell_Activated"            = "#8B0000",
  "Aerocyte_Endothelial"           = "#2CA02C",
  "Bcell_Naive"                    = "#637939",
  "Plasma_Cell_MZB1"               = "#0000FF",
  "Fibroblast_Matrix"              = "#BCBD22",
  "Vascular_Endothelial"           = "#3182BD",
  "Pericyte_SMC"                   = "#17BECF",
  "Mast_Cell"                      = "#E6550D",
  "Alveolar_Macrophage"            = "#66A61E",
  "Club_Ciliated_Transitional"     = "#FF00FF",
  "Epithelial_Progenitor"          = "#E377C2",
  "Goblet_MUC5B"                   = "#393B79",
  "AT2_Transitional"               = "#E6AB02"
)

# --------------------------------
# Convert cluster IDs to harmonized cell types
# --------------------------------

segmap_data <- segmap_data %>%
  mutate(
    cluster = as.character(cluster),
    
    cell_type = factor(
      cluster_annotation_segmap_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
# --------------------------------

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
segmap_data_spatial <- readr::read_csv("XeniumData_segmap/segmapOutput_4/cellpose_pipeline_outputs_1_iter15/data/spatial_clusters_final.csv", show_col_types = FALSE)
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