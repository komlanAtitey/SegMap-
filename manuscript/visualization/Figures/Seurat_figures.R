library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
# ==============================================================================
#
#
# DATA 1
#
#
# ==============================================================================
#===================================================
# UMAP embedding 1
#===================================================
load("XeniumData_seurat/new_xenium_data1.rdata")

cluster_annotation_seurat_harmonized <- c(
  "1" = "Oligodendrocyte_Reactive_SPP1",
  "2" = "Oligodendrocyte_Stressed_CRYAB",
  "3" = "Excitatory_IT_MEF2C",
  "4" = "Endothelial_FLT1",
  "5" = "Astrocyte_Reactive_GFAP",
  "6" = "Microglia_DAM_SPP1",
  "7" = "Interneuron_LAMP5_KIT",
  "8" = "OPC_PDGFRA"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Oligodendrocyte_Reactive_SPP1" = "#0000FF",
  "Oligodendrocyte_Stressed_CRYAB"= "#FF0000",
  "Excitatory_IT_MEF2C"           = "#31A354",
  "Endothelial_FLT1"              = "#E377C2",
  "Astrocyte_Reactive_GFAP"       = "#00FF7F",
  "Microglia_DAM_SPP1"            = "#8C564B",
  "Interneuron_LAMP5_KIT"         = "#E6550D",
  "OPC_PDGFRA"                    = "#17BECF"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

seurat_data <- as.data.frame(
  Embeddings(new_xenium_data1, reduction = "umap")
)

colnames(seurat_data) <- c("UMAP_1", "UMAP_2")

seurat_data$cluster <- as.character(Idents(new_xenium_data1))

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

seurat_data <- seurat_data %>%
  mutate(
    cell_type = factor(
      cluster_annotation_seurat_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Seurat Cell Types",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 15)

p_seurat


#===================================================
# spatial visualization 1
#===================================================

# --------------------------------
# Extract spatial coordinates
# --------------------------------

coords <- GetTissueCoordinates(new_xenium_data1)

coords$cluster <- as.character(
  Idents(new_xenium_data1)
)

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

coords$cell_type <- factor(
  cluster_annotation_seurat_harmonized[coords$cluster],
  levels = names(celltype_colors)
)

# --------------------------------
# Spatial Plot
# --------------------------------

p_seurat_spatial <- ggplot(
  coords,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  scale_y_reverse() +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X Coordinate",
    y = "Y Coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_seurat_spatial

# Without legend
p_seurat_spatial + NoLegend()

# ==============================================================================
#
#
# DATA 2
#
#
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)

load("XeniumData_seurat/new_xenium_data2.rdata")
# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_seurat_harmonized <- c(
  "1"  = "LAMP5_CUX2_Interneuron",
  "2"  = "Oligodendrocyte_Mature",
  "3"  = "Astrocyte_Reactive_APOE",
  "4"  = "Endothelial_FLT1",
  "5"  = "PVALB_TAC1_Interneuron",
  "6"  = "Excitatory_RORB_HS3ST2",
  "7"  = "Microglia_DAM_ITGAX",
  "8"  = "Pericyte_ABCC9",
  "9"  = "VIP_RELN_Interneuron",
  "10" = "OPC_PDGFRA",
  "11" = "KIT_RELN_LAMP5_Interneuron"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "LAMP5_CUX2_Interneuron"      = "#8C564B",
  "Oligodendrocyte_Mature"      = "#2CA02C",
  "Astrocyte_Reactive_APOE"     = "#3182BD",
  "Endothelial_FLT1"            = "#0000FF",
  "PVALB_TAC1_Interneuron"      = "#66A61E",
  "Excitatory_RORB_HS3ST2"      = "#551A8B",
  "Microglia_DAM_ITGAX"         = "#CD7054",
  "Pericyte_ABCC9"              = "#17BECF",
  "VIP_RELN_Interneuron"        = "#A6761D",
  "OPC_PDGFRA"                  = "#31A354",
  "KIT_RELN_LAMP5_Interneuron"  = "#FF1493"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

seurat_data <- as.data.frame(
  Embeddings(new_xenium_data2, reduction = "umap")
)

colnames(seurat_data) <- c("UMAP_1", "UMAP_2")

seurat_data$cluster <- as.character(
  Idents(new_xenium_data2)
)

# --------------------------------
# Check for unmapped clusters
# --------------------------------

missing_clusters <- setdiff(
  unique(seurat_data$cluster),
  names(cluster_annotation_seurat_harmonized)
)

if(length(missing_clusters) > 0){
  message(
    "Clusters missing annotation: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

seurat_data <- seurat_data %>%
  mutate(
    cell_type = cluster_annotation_seurat_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
    size = 0.30,
    alpha = 0.80
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Seurat Cell Types",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
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

p_seurat
# p_seurat + NoLegend()

#===================================================
# spatial visualization 2
#===================================================

# --------------------------------
# Extract spatial coordinates
# --------------------------------

coords <- GetTissueCoordinates(new_xenium_data2)

coords$cluster <- as.character(
  Idents(new_xenium_data2)
)

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

coords$cell_type <- factor(
  cluster_annotation_seurat_harmonized[coords$cluster],
  levels = names(celltype_colors)
)

# --------------------------------
# Spatial Plot
# --------------------------------

p_seurat_spatial <- ggplot(
  coords,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  scale_y_reverse() +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X Coordinate",
    y = "Y Coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_seurat_spatial

# Without legend
p_seurat_spatial + NoLegend()

# ==============================================================================
#
#
# DATA 3
#
#
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)

load("XeniumData_seurat/new_xenium_data3.rdata")
# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_seurat_harmonized <- c(
  "1"  = "Interneuron_LAMP5_CUX2",
  "2"  = "Oligodendrocyte_Mature",
  "3"  = "Astrocyte_Reactive",
  "4"  = "Endothelial_FLT1",
  "5"  = "Interneuron_PVALB_LHX6",
  "6"  = "Excitatory_TSHZ2_HS3ST2",
  "7"  = "Microglia_DAM",
  "8"  = "Pericyte_ABCC9",
  "9"  = "Interneuron_VIP_PROX1",
  "10" = "OPC_PDGFRA",
  "11" = "Interneuron_KIT_RELN"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Interneuron_LAMP5_CUX2"    = "#637939",
  "Oligodendrocyte_Mature"    = "#2CA02C",
  "Astrocyte_Reactive"        = "#3182BD",
  "Endothelial_FLT1"          = "#0000FF",
  "Interneuron_PVALB_LHX6"    = "#756BB1",
  "Excitatory_TSHZ2_HS3ST2"   = "#17BECF",
  "Microglia_DAM"             = "#FFE4E1",
  "Pericyte_ABCC9"            = "#00FF7F",
  "Interneuron_VIP_PROX1"     = "#8B0000",
  "OPC_PDGFRA"                = "#31A354",
  "Interneuron_KIT_RELN"      = "#CD7054"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

seurat_data <- as.data.frame(
  Embeddings(new_xenium_data3, reduction = "umap")
)

colnames(seurat_data) <- c("UMAP_1", "UMAP_2")

seurat_data$cluster <- as.character(
  Idents(new_xenium_data3)
)

# --------------------------------
# Check for unmapped clusters
# --------------------------------

missing_clusters <- setdiff(
  unique(seurat_data$cluster),
  names(cluster_annotation_seurat_harmonized)
)

if(length(missing_clusters) > 0){
  message(
    "Clusters missing annotation: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

seurat_data <- seurat_data %>%
  mutate(
    cell_type = cluster_annotation_seurat_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
    size = 0.30,
    alpha = 0.80
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Seurat Cell Types",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
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

p_seurat
# p_seurat + NoLegend()

#===================================================
# spatial visualization 2
#===================================================

# --------------------------------
# Extract spatial coordinates
# --------------------------------

coords <- GetTissueCoordinates(new_xenium_data3)

coords$cluster <- as.character(
  Idents(new_xenium_data3)
)

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

coords$cell_type <- factor(
  cluster_annotation_seurat_harmonized[coords$cluster],
  levels = names(celltype_colors)
)

# --------------------------------
# Spatial Plot
# --------------------------------

p_seurat_spatial <- ggplot(
  coords,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  scale_y_reverse() +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X Coordinate",
    y = "Y Coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_seurat_spatial

# Without legend
p_seurat_spatial + NoLegend()

# ==============================================================================
#
#
# DATA 4
#
#
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)

load("XeniumData_seurat/new_xenium_data4.rdata")
# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_seurat_harmonized <- c(
  "1"  = "CD4_Tcell_Activated",
  "2"  = "Epithelial_Progenitor",
  "3"  = "Aerocyte_Endothelial",
  "4"  = "Macrophage_FCGR3A",
  "5"  = "Fibroblast_Matrix",
  "6"  = "Pericyte_SMC",
  "7"  = "Goblet_MUC5B",
  "8"  = "Bcell_Naive",
  "9"  = "AT2_Alveolar",
  "10" = "Plasma_Cell_MZB1",
  "11" = "AT1_Alveolar",
  "12" = "Mast_Cell",
  "13" = "Lymphatic_Endothelial"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "CD4_Tcell_Activated"          = "#8B0000",
  "Epithelial_Progenitor"        = "#E377C2",
  "Aerocyte_Endothelial"         = "#2CA02C",
  "Macrophage_FCGR3A"            = "#00FFFF",
  "Fibroblast_Matrix"            = "#BCBD22",
  "Pericyte_SMC"                 = "#17BECF",
  "Goblet_MUC5B"                 = "#393B79",
  "Bcell_Naive"                  = "#637939",
  "AT2_Alveolar"                 = "#551A8B",
  "Plasma_Cell_MZB1"             = "#0000FF",
  "AT1_Alveolar"                 = "#1B9E77",
  "Mast_Cell"                    = "#E6550D",
  "Lymphatic_Endothelial"        = "#FFE4E1"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

seurat_data <- as.data.frame(
  Embeddings(new_xenium_data4, reduction = "umap")
)

colnames(seurat_data) <- c("UMAP_1", "UMAP_2")

seurat_data$cluster <- as.character(
  Idents(new_xenium_data4)
)

# --------------------------------
# Check for unmapped clusters
# --------------------------------

missing_clusters <- setdiff(
  unique(seurat_data$cluster),
  names(cluster_annotation_seurat_harmonized)
)

if(length(missing_clusters) > 0){
  message(
    "Clusters missing annotation: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

seurat_data <- seurat_data %>%
  mutate(
    cell_type = cluster_annotation_seurat_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
    size = 0.30,
    alpha = 0.80
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Seurat Cell Types",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
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

p_seurat
# p_seurat + NoLegend()

#===================================================
# spatial visualization 2
#===================================================

# --------------------------------
# Extract spatial coordinates
# --------------------------------

coords <- GetTissueCoordinates(new_xenium_data4)

coords$cluster <- as.character(
  Idents(new_xenium_data4)
)

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

coords$cell_type <- factor(
  cluster_annotation_seurat_harmonized[coords$cluster],
  levels = names(celltype_colors)
)

# --------------------------------
# Spatial Plot
# --------------------------------

p_seurat_spatial <- ggplot(
  coords,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  scale_y_reverse() +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Cell Types",
    x = "X Coordinate",
    y = "Y Coordinate"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.key = element_blank(),
    legend.position = "right",
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  )

p_seurat_spatial

# Without legend
p_seurat_spatial + NoLegend()

