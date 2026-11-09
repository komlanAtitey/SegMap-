

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
load("XeniumData_baysor/baysor_data_1.rdata")

cluster_annotation_baysor_harmonized <- c(
  "0" = "Astrocyte_Reactive_GFAP",
  "1" = "Oligodendrocyte_Reactive_SPP1",
  "2" = "Excitatory_IT_MEF2C",
  "3" = "Endothelial_FLT1",
  "4" = "Oligodendrocyte_Mature",
  "5" = "OPC_PDGFRA",
  "6" = "Unclassified_Stress_Response_a",
  "7" = "Unclassified_Stress_Response_b"
)

# --------------------------------
# Cell-type colors
# --------------------------------
celltype_colors <- c(
  "Astrocyte_Reactive_GFAP"         = "#00FF7F",
  "Oligodendrocyte_Reactive_SPP1"   = "#0000FF",
  "Excitatory_IT_MEF2C"             = "#31A354",
  "Endothelial_FLT1"                = "#E377C2",
  "Oligodendrocyte_Mature"          = "#1F77B4",
  "OPC_PDGFRA"                      = "#17BECF",
  "Unclassified_Stress_Response_a"  = "#1B9E77",
  "Unclassified_Stress_Response_b"  = "#637939"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

baysor_data <- as.data.frame(
  Embeddings(baysor_data_1, reduction = "umap")
)

colnames(baysor_data) <- c("UMAP_1", "UMAP_2")

baysor_data$cluster <- as.character(Idents(baysor_data_1))

# --------------------------------
# Convert clusters to harmonized cell types
# --------------------------------

baysor_data <- baysor_data %>%
  mutate(
    cell_type = factor(
      cluster_annotation_baysor_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
# --------------------------------

p_segmap <- ggplot(
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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of SegMap Cell Types",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 15)

p_segmap


#===================================================
# spatial visualization 1
#===================================================

# --------------------------------
# Create plotting dataframe
# --------------------------------
baysor_spatial <- baysor_data_1@meta.data

baysor_spatial$cluster <- as.character(
  baysor_spatial$seurat_clusters
)

baysor_spatial <- baysor_spatial %>%
  mutate(
    cell_type = factor(
      cluster_annotation_baysor_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# Check for missing annotations
# --------------------------------
table(baysor_spatial$cell_type, useNA = "ifany")

setdiff(
  unique(baysor_spatial$cluster),
  names(cluster_annotation_baysor_harmonized)
)

# --------------------------------
# Spatial plot
# --------------------------------
p_baysor_spatial <- ggplot(
  baysor_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.20,
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
    title = "Baysor Spatial Cell-Type Map",
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

p_baysor_spatial

# ==============================================================================
#
#
# DATA 2
#
#
# ==============================================================================
#===================================================
# UMAP embedding 1
#===================================================
load("XeniumData_baysor/baysor_data_2.rdata")

library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_baysor_harmonized <- c(
  "0"  = "Astrocyte_Reactive_APOE",
  "1"  = "Excitatory_TSHZ2_CRYM",
  "2"  = "Oligodendrocyte_Mature",
  "3"  = "Astrocyte_Homeostatic_FGFR3",
  "4"  = "Endothelial_FLT1",
  "5"  = "SST_PVALB_LHX6_Interneuron",
  "6"  = "VIP_PROX1_RELN_Interneuron",
  "7"  = "Microglia_DAM_TREM2",
  "8"  = "Excitatory_LMO4_NRGN",
  "9"  = "OPC_PDGFRA_GPNMB",
  "10" = "Excitatory_IDH2_NRGN"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Astrocyte_Reactive_APOE"       = "#3182BD",
  "Excitatory_TSHZ2_CRYM"         = "#FF7F0E",
  "Oligodendrocyte_Mature"        = "#2CA02C",
  "Astrocyte_Homeostatic_FGFR3"   = "#00FF7F",
  "Endothelial_FLT1"              = "#0000FF",
  "SST_PVALB_LHX6_Interneuron"    = "#7B4173",
  "VIP_PROX1_RELN_Interneuron"    = "#551A8B",
  "Microglia_DAM_TREM2"           = "#BCBD22",
  "Excitatory_LMO4_NRGN"          = "#66A61E",
  "OPC_PDGFRA_GPNMB"              = "#FF0000",
  "Excitatory_IDH2_NRGN"          = "#E377C2"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

baysor_data <- as.data.frame(
  Embeddings(baysor_data_2, reduction = "umap")
)

colnames(baysor_data) <- c("UMAP_1", "UMAP_2")

baysor_data$cluster <- as.character(
  Idents(baysor_data_2)
)

# --------------------------------
# Check for missing annotations
# --------------------------------

missing_clusters <- setdiff(
  unique(baysor_data$cluster),
  names(cluster_annotation_baysor_harmonized)
)

if(length(missing_clusters) > 0){
  warning(
    "Missing annotations for clusters: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

baysor_data <- baysor_data %>%
  mutate(
    cell_type = cluster_annotation_baysor_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Baysor Cell Types",
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
    legend.text = element_text(
      size = 10
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

p_baysor
# p_baysor + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

# --------------------------------
# Create plotting dataframe
# --------------------------------
baysor_spatial <- baysor_data_2@meta.data

baysor_spatial$cluster <- as.character(
  baysor_spatial$seurat_clusters
)

baysor_spatial <- baysor_spatial %>%
  mutate(
    cell_type = factor(
      cluster_annotation_baysor_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# Check for missing annotations
# --------------------------------
table(baysor_spatial$cell_type, useNA = "ifany")

setdiff(
  unique(baysor_spatial$cluster),
  names(cluster_annotation_baysor_harmonized)
)

# --------------------------------
# Spatial plot
# --------------------------------
p_baysor_spatial <- ggplot(
  baysor_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.20,
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
    title = "Baysor Spatial Cell-Type Map",
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

p_baysor_spatial

# ==============================================================================
#
#
# DATA 3
#
#
# ==============================================================================
#===================================================
# UMAP embedding 1
#===================================================
load("XeniumData_baysor/baysor_data_3.rdata")

library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_baysor_harmonized <- c(
  "0"  = "Excitatory_CUX2_RORB",
  "1"  = "Excitatory_CRYM_EPHA4",
  "2"  = "Oligodendrocyte_Mature",
  "3"  = "Astrocyte_Reactive",
  "4"  = "Endothelial_FLT1",
  "5"  = "Interneuron_SST_PVALB",
  "6"  = "Interneuron_VIP_PROX1",
  "7"  = "Microglia_DAM",
  "8"  = "Excitatory_NRGN_LMO4",
  "9"  = "OPC_PDGFRA",
  "10" = "Excitatory_LMO4_NPNT"
)

# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Excitatory_CUX2_RORB"      = "#BCBD22",
  "Excitatory_CRYM_EPHA4"     = "#A6761D",
  "Oligodendrocyte_Mature"    = "#2CA02C",
  "Astrocyte_Reactive"        = "#3182BD",
  "Endothelial_FLT1"          = "#0000FF",
  "Interneuron_SST_PVALB"     = "#66A61E",
  "Interneuron_VIP_PROX1"     = "#8B0000",
  "Microglia_DAM"             = "#FFE4E1",
  "Excitatory_NRGN_LMO4"      = "#A6761D",
  "OPC_PDGFRA"                = "#31A354",
  "Excitatory_LMO4_NPNT"      = "#393B79"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

baysor_data <- as.data.frame(
  Embeddings(baysor_data_3, reduction = "umap")
)

colnames(baysor_data) <- c("UMAP_1", "UMAP_2")

baysor_data$cluster <- as.character(
  Idents(baysor_data_3)
)

# --------------------------------
# Check for missing annotations
# --------------------------------

missing_clusters <- setdiff(
  unique(baysor_data$cluster),
  names(cluster_annotation_baysor_harmonized)
)

if(length(missing_clusters) > 0){
  warning(
    "Missing annotations for clusters: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

baysor_data <- baysor_data %>%
  mutate(
    cell_type = cluster_annotation_baysor_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Baysor Cell Types",
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
    legend.text = element_text(
      size = 10
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

p_baysor
# p_baysor + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

# --------------------------------
# Create plotting dataframe
# --------------------------------
baysor_spatial <- baysor_data_3@meta.data

baysor_spatial$cluster <- as.character(
  baysor_spatial$seurat_clusters
)

baysor_spatial <- baysor_spatial %>%
  mutate(
    cell_type = factor(
      cluster_annotation_baysor_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# Check for missing annotations
# --------------------------------
table(baysor_spatial$cell_type, useNA = "ifany")

setdiff(
  unique(baysor_spatial$cluster),
  names(cluster_annotation_baysor_harmonized)
)

# --------------------------------
# Spatial plot
# --------------------------------
p_baysor_spatial <- ggplot(
  baysor_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.20,
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
    title = "Baysor Spatial Cell-Type Map",
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

p_baysor_spatial

# ==============================================================================
#
#
# DATA 4
#
#
# ==============================================================================
#===================================================
# UMAP embedding 1
#===================================================
load("XeniumData_baysor/baysor_data_4.rdata")

library(Seurat)
library(dplyr)
library(ggplot2)

# --------------------------------
# Harmonized annotations
# --------------------------------

cluster_annotation_baysor_harmonized <- c(
  "0"  = "Secretory_Epithelial_ACE2",
  "1"  = "Macrophage_MARCO",
  "2"  = "CD8_Tcell_Cytotoxic",
  "3"  = "Fibroblast_PDGFRA",
  "4"  = "Aerocyte_Endothelial",
  "5"  = "Vascular_Endothelial",
  "6"  = "AT1_Alveolar",
  "7"  = "Pericyte_MYH11",
  "8"  = "Plasma_Cell_MZB1",
  "9"  = "Ciliated_Epithelial",
  "10" = "Epithelial_Progenitor",
  "11" = "Mast_Cell",
  "12" = "Lymphatic_Endothelial"
)


# --------------------------------
# Cell-type colors
# --------------------------------

celltype_colors <- c(
  "Secretory_Epithelial_ACE2"      = "#17BECF",
  "Macrophage_MARCO"               = "#FF0000",
  "CD8_Tcell_Cytotoxic"            = "#CD7054",
  "Fibroblast_PDGFRA"              = "#7B4173",
  "Aerocyte_Endothelial"           = "#2CA02C",
  "Vascular_Endothelial"           = "#3182BD",
  "AT1_Alveolar"                   = "#1B9E77",
  "Pericyte_MYH11"                 = "#756BB1",
  "Plasma_Cell_MZB1"               = "#0000FF",
  "Ciliated_Epithelial"            = "#31A354",
  "Epithelial_Progenitor"          = "#E377C2",
  "Mast_Cell"                      = "#E6550D",
  "Lymphatic_Endothelial"          = "#FFE4E1"
)

# --------------------------------
# Extract UMAP coordinates
# --------------------------------

baysor_data <- as.data.frame(
  Embeddings(baysor_data_4, reduction = "umap")
)

colnames(baysor_data) <- c("UMAP_1", "UMAP_2")

baysor_data$cluster <- as.character(
  Idents(baysor_data_4)
)

# --------------------------------
# Check for missing annotations
# --------------------------------

missing_clusters <- setdiff(
  unique(baysor_data$cluster),
  names(cluster_annotation_baysor_harmonized)
)

if(length(missing_clusters) > 0){
  warning(
    "Missing annotations for clusters: ",
    paste(missing_clusters, collapse = ", ")
  )
}

# --------------------------------
# Convert clusters to cell types
# --------------------------------

baysor_data <- baysor_data %>%
  mutate(
    cell_type = cluster_annotation_baysor_harmonized[cluster],
    cell_type = factor(
      cell_type,
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# UMAP Plot
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
  scale_color_manual(
    values = celltype_colors,
    drop = FALSE,
    na.value = "black",
    name = "Cell Type"
  ) +
  coord_equal() +
  labs(
    title = "UMAP of Baysor Cell Types",
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
    legend.text = element_text(
      size = 10
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

p_baysor
# p_baysor + NoLegend()

#===================================================
# spatial visualization 1
#===================================================

# --------------------------------
# Create plotting dataframe
# --------------------------------
baysor_spatial <- baysor_data_4@meta.data

baysor_spatial$cluster <- as.character(
  baysor_spatial$seurat_clusters
)

baysor_spatial <- baysor_spatial %>%
  mutate(
    cell_type = factor(
      cluster_annotation_baysor_harmonized[cluster],
      levels = names(celltype_colors)
    )
  )

# --------------------------------
# Check for missing annotations
# --------------------------------
table(baysor_spatial$cell_type, useNA = "ifany")

setdiff(
  unique(baysor_spatial$cluster),
  names(cluster_annotation_baysor_harmonized)
)

# --------------------------------
# Spatial plot
# --------------------------------
p_baysor_spatial <- ggplot(
  baysor_spatial,
  aes(
    x = x,
    y = y,
    color = cell_type
  )
) +
  geom_point(
    size = 0.20,
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
    title = "Baysor Spatial Cell-Type Map",
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

p_baysor_spatial