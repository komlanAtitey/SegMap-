library(dplyr)
library(ggplot2)
library(readr)
library(grid)

# ==============================================================================
#
#
# DATA 1
#
#
# ==============================================================================

#===================================================
# umap embedding 1
#===================================================
#---------------------------
# Load data
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_1/umap_with_clusters.csv",
  show_col_types = FALSE
)

# Check column names
print(colnames(scanpy_data))

#---------------------------
# Harmonized annotations
#---------------------------

cluster_annotation_scanpy_harmonized <- c(
  "0" = "Oligodendrocyte_Reactive_SPP1",
  "1" = "Oligodendrocyte_MOBP_Stressed",
  "2" = "Endothelial_FLT1",
  "3" = "Astrocyte_Homeostatic",
  "4" = "Microglia_DAM_SPP1",
  "5" = "Excitatory_IT_MEF2C",
  "6" = "Immune_PTPRC",
  "7" = "Astrocyte_Reactive_AQP4"
)

#---------------------------
# Cell-type colors
#---------------------------

celltype_colors <- c(
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Oligodendrocyte_MOBP_Stressed" = "#FF00FF",
  "Endothelial_FLT1"              = "#E377C2",
  "Astrocyte_Homeostatic"         = "#D62728",
  "Microglia_DAM_SPP1"            = "#8C564B",
  "Excitatory_IT_MEF2C"           = "#31A354",
  "Immune_PTPRC"                  = "#637939",
  "Astrocyte_Reactive_AQP4"       = "#66A61E"
)

#---------------------------
# Rename columns if needed
#---------------------------

if ("leiden" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, cluster = leiden)
}

if ("UMAP1" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_1 = UMAP1)
}

if ("UMAP2" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_2 = UMAP2)
}

#---------------------------
# Convert clusters to cell types
#---------------------------

scanpy_data <- scanpy_data %>%
  mutate(
    cluster = as.character(cluster),
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# UMAP Plot
#---------------------------

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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  labs(
    title = "UMAP of Scanpy Cell Types",
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
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

p_scanpy + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

#---------------------------
# Load Scanpy UMAP clusters
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_1/umap_with_clusters.csv",
  show_col_types = FALSE
)

#---------------------------
# Load spatial coordinates
#---------------------------

scanpy_data_spatial <- read_csv(
  "XeniumData_scanpy/scanpy_data_1/spatial_coords.csv",
  show_col_types = FALSE
)

# Inspect columns
print(colnames(scanpy_data_spatial))

#---------------------------
# Add cluster labels
#---------------------------

scanpy_data_spatial$cluster <- as.character(scanpy_data$cluster)

#---------------------------
# Convert clusters to harmonized cell types
#---------------------------

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# Spatial Plot
#---------------------------
p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "Scanpy Spatial Cell-Type Map",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
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
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_scanpy_spatial
p_scanpy_spatial + NoLegend()

# Remove legend if desired
# p_scanpy_spatial + NoLegend()
# ==============================================================================
#
#
# DATA 2
#
#
# ==============================================================================

#===================================================
# umap embedding 1
#===================================================
#---------------------------
# Load data
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_2/umap_with_clusters.csv",
  show_col_types = FALSE
)

# Check column names
print(colnames(scanpy_data))

#---------------------------
# Harmonized annotations
#---------------------------
cluster_annotation_scanpy_harmonized <- c(
  "0"  = "Astrocyte_Reactive_APOE",
  "1"  = "PVALB_TAC1_Interneuron",
  "2"  = "Excitatory_HTR2A_CCK",
  "3"  = "Excitatory_LAMP5_C1QL3",
  "4"  = "Excitatory_RORB_HS3ST2",
  "5"  = "Oligodendrocyte_Mature",
  "6"  = "Endothelial_FLT1",
  "7"  = "OPC_PDGFRA",
  "8"  = "Microglia_DAM_ITGAX",
  "9"  = "VIP_RELN_Interneuron",
  "10" = "KIT_RELN_LAMP5_Interneuron"
)

#---------------------------
# Cell-type colors
#---------------------------

celltype_colors <- c(
  "Astrocyte_Reactive_APOE"     = "#3182BD",
  "PVALB_TAC1_Interneuron"      = "#66A61E",
  "Excitatory_HTR2A_CCK"        = "#E6AB02",
  "Excitatory_LAMP5_C1QL3"      = "#E6550D",
  "Excitatory_RORB_HS3ST2"      = "#551A8B",
  "Oligodendrocyte_Mature"      = "#2CA02C",
  "Endothelial_FLT1"            = "#0000FF",
  "OPC_PDGFRA"                  = "#31A354",
  "Microglia_DAM_ITGAX"         = "#CD7054",
  "VIP_RELN_Interneuron"        = "#A6761D",
  "KIT_RELN_LAMP5_Interneuron"  = "#FF1493"
)
#---------------------------
# Rename columns if needed
#---------------------------

if ("leiden" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, cluster = leiden)
}

if ("UMAP1" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_1 = UMAP1)
}

if ("UMAP2" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_2 = UMAP2)
}

#---------------------------
# Convert clusters to cell types
#---------------------------

scanpy_data <- scanpy_data %>%
  mutate(
    cluster = as.character(cluster),
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# UMAP Plot
#---------------------------

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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  labs(
    title = "UMAP of Scanpy Cell Types",
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
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

p_scanpy + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

#---------------------------
# Load Scanpy UMAP clusters
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_2/umap_with_clusters.csv",
  show_col_types = FALSE
)

#---------------------------
# Load spatial coordinates
#---------------------------

scanpy_data_spatial <- read_csv(
  "XeniumData_scanpy/scanpy_data_2/spatial_coords.csv",
  show_col_types = FALSE
)

# Inspect columns
print(colnames(scanpy_data_spatial))

#---------------------------
# Add cluster labels
#---------------------------

scanpy_data_spatial$cluster <- as.character(scanpy_data$cluster)

#---------------------------
# Convert clusters to harmonized cell types
#---------------------------

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# Spatial Plot
#---------------------------
p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "Scanpy Spatial Cell-Type Map",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
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
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_scanpy_spatial
p_scanpy_spatial + NoLegend()

# ==============================================================================
#
#
# DATA 3
#
#
# ==============================================================================

#===================================================
# umap embedding 1
#===================================================
#---------------------------
# Load data
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_3/umap_with_clusters.csv",
  show_col_types = FALSE
)

# Check column names
print(colnames(scanpy_data))

#---------------------------
# Harmonized annotations
#---------------------------
cluster_annotation_scanpy_harmonized <- c(
  "0"  = "Astrocyte_Reactive",
  "1"  = "Interneuron_PVALB_LHX6",
  "2"  = "Excitatory_CCK_HTR2A",
  "3"  = "Excitatory_LAMP5_C1QL3",
  "4"  = "Excitatory_RORB_HS3ST2",
  "5"  = "Oligodendrocyte_Mature",
  "6"  = "Endothelial_FLT1",
  "7"  = "OPC_PDGFRA",
  "8"  = "Microglia_DAM",
  "9"  = "Interneuron_VIP_PROX1",
  "10" = "Interneuron_KIT_LAMP5"
)

#---------------------------
# Cell-type colors
#---------------------------
celltype_colors <- c(
  "Astrocyte_Reactive"       = "#3182BD",
  "Interneuron_PVALB_LHX6"   = "#756BB1",
  "Excitatory_CCK_HTR2A"     = "#E6AB02",
  "Excitatory_LAMP5_C1QL3"   = "#E6550D",
  "Excitatory_RORB_HS3ST2"   = "#551A8B",
  "Oligodendrocyte_Mature"   = "#2CA02C",
  "Endothelial_FLT1"         = "#0000FF",
  "OPC_PDGFRA"               = "#31A354",
  "Microglia_DAM"            = "#FFE4E1",
  "Interneuron_VIP_PROX1"    = "#8B0000",
  "Interneuron_KIT_LAMP5"    = "#00FFFF"
)

#---------------------------
# Rename columns if needed
#---------------------------

if ("leiden" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, cluster = leiden)
}

if ("UMAP1" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_1 = UMAP1)
}

if ("UMAP2" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_2 = UMAP2)
}

#---------------------------
# Convert clusters to cell types
#---------------------------

scanpy_data <- scanpy_data %>%
  mutate(
    cluster = as.character(cluster),
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# UMAP Plot
#---------------------------

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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  labs(
    title = "UMAP of Scanpy Cell Types",
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
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

p_scanpy + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

#---------------------------
# Load Scanpy UMAP clusters
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_3/umap_with_clusters.csv",
  show_col_types = FALSE
)

#---------------------------
# Load spatial coordinates
#---------------------------

scanpy_data_spatial <- read_csv(
  "XeniumData_scanpy/scanpy_data_3/spatial_coords.csv",
  show_col_types = FALSE
)

# Inspect columns
print(colnames(scanpy_data_spatial))

#---------------------------
# Add cluster labels
#---------------------------

scanpy_data_spatial$cluster <- as.character(scanpy_data$cluster)

#---------------------------
# Convert clusters to harmonized cell types
#---------------------------

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# Spatial Plot
#---------------------------
p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "Scanpy Spatial Cell-Type Map",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
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
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_scanpy_spatial
p_scanpy_spatial + NoLegend()

# ==============================================================================
#
#
# DATA 4
#
#
# ==============================================================================

#===================================================
# umap embedding 1
#===================================================
#---------------------------
# Load data
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_4/umap_with_clusters.csv",
  show_col_types = FALSE
)

# Check column names
print(colnames(scanpy_data))

#---------------------------
# Harmonized annotations
#---------------------------

cluster_annotation_scanpy_harmonized <- c(
  "0"  = "Fibroblast_PDGFRA",
  "1"  = "Pericyte_MYH11",
  "2"  = "CD4_Tcell_Activated",
  "3"  = "Vascular_Endothelial",
  "4"  = "Macrophage_FCGR3A",
  "5"  = "Plasma_Cell_MZB1",
  "6"  = "Epithelial_Progenitor",
  "7"  = "Bcell_Naive",
  "8"  = "AT1_Alveolar",
  "9"  = "Aerocyte_Endothelial",
  "10" = "AT2_Alveolar",
  "11" = "Goblet_MUC5B",
  "12" = "Mast_Cell"
)

#---------------------------
# Cell-type colors
#---------------------------

celltype_colors <- c(
  "Fibroblast_PDGFRA"              = "#7B4173",
  "Pericyte_MYH11"                 = "#756BB1",
  "CD4_Tcell_Activated"            = "#8B0000",
  "Vascular_Endothelial"           = "#3182BD",
  "Macrophage_FCGR3A"              = "#00FFFF",
  "Plasma_Cell_MZB1"               = "#0000FF",
  "Epithelial_Progenitor"          = "#E377C2",
  "Bcell_Naive"                    = "#637939",
  "AT1_Alveolar"                   = "#1B9E77",
  "Aerocyte_Endothelial"           = "#2CA02C",
  "AT2_Alveolar"                   = "#551A8B",
  "Goblet_MUC5B"                   = "#393B79",
  "Mast_Cell"                      = "#E6550D"
)

#---------------------------
# Rename columns if needed
#---------------------------

if ("leiden" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, cluster = leiden)
}

if ("UMAP1" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_1 = UMAP1)
}

if ("UMAP2" %in% colnames(scanpy_data)) {
  scanpy_data <- rename(scanpy_data, UMAP_2 = UMAP2)
}

#---------------------------
# Convert clusters to cell types
#---------------------------

scanpy_data <- scanpy_data %>%
  mutate(
    cluster = as.character(cluster),
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# UMAP Plot
#---------------------------

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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  labs(
    title = "UMAP of Scanpy Cell Types",
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
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

p_scanpy #+ NoLegend()

#===================================================
# spatial visualization 1
#===================================================

#---------------------------
# Load Scanpy UMAP clusters
#---------------------------

scanpy_data <- read_csv(
  "XeniumData_scanpy/scanpy_data_4/umap_with_clusters.csv",
  show_col_types = FALSE
)

#---------------------------
# Load spatial coordinates
#---------------------------

scanpy_data_spatial <- read_csv(
  "XeniumData_scanpy/scanpy_data_4/spatial_coords.csv",
  show_col_types = FALSE
)

# Inspect columns
print(colnames(scanpy_data_spatial))

#---------------------------
# Add cluster labels
#---------------------------

scanpy_data_spatial$cluster <- as.character(scanpy_data$cluster)

#---------------------------
# Convert clusters to harmonized cell types
#---------------------------

scanpy_data_spatial <- scanpy_data_spatial %>%
  mutate(
    cell_type = cluster_annotation_scanpy_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

#---------------------------
# Spatial Plot
#---------------------------
p_scanpy_spatial <- ggplot(
  scanpy_data_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.25,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title = "Scanpy Spatial Cell-Type Map",
    x = "X coordinate",
    y = "Y coordinate"
  ) +
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
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_scanpy_spatial
p_scanpy_spatial + NoLegend()


