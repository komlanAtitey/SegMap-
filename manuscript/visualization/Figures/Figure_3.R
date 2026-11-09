setwd("/Users/kxa643/Documents/xenium")
getwd()

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tidyverse)
library(viridis)
library(tibble)
library(readr)
library(purrr)
library(knitr)
library(Matrix)
library(Seurat)
library(spdep)

# ==============================================================================
#      Figure 3a: Mean -log10(FDR),  Figure 2b: Significant Marker Genes
# ==============================================================================

# ----------------------
# Add dataset labels
# ----------------------
data_DE_per_cluster_1 <- read.csv("segmap_data_Figure_3/data_DE_per_cluster_1.csv", header = TRUE, stringsAsFactors = FALSE)
data_DE_per_cluster_2 <- read.csv("segmap_data_Figure_3/data_DE_per_cluster_2.csv", header = TRUE, stringsAsFactors = FALSE)
data_DE_per_cluster_3 <- read.csv("segmap_data_Figure_3/data_DE_per_cluster_3.csv", header = TRUE, stringsAsFactors = FALSE)
data_DE_per_cluster_4 <- read.csv("segmap_data_Figure_3/data_DE_per_cluster_4.csv", header = TRUE, stringsAsFactors = FALSE)

data_DE_per_cluster_1$Dataset <- "Dataset 1"
data_DE_per_cluster_2$Dataset <- "Dataset 2"
data_DE_per_cluster_3$Dataset <- "Dataset 3"
data_DE_per_cluster_4$Dataset <- "Dataset 4"

# ----------------------
# Merge all datasets
# ----------------------
benchmark <- bind_rows(
  data_DE_per_cluster_1,
  data_DE_per_cluster_2,
  data_DE_per_cluster_3,
  data_DE_per_cluster_4
)

benchmark$method <- factor(
  benchmark$method,
  levels=c("Segmap","Scanpy","Seurat","Baysor")
)

# ----------------------
# Compute one value per dataset and method
# ----------------------
summary.dataset <- benchmark %>%
  group_by(Dataset, method) %>%
  summarise(
    
    median_log2FC =
      median(median_log2FC, na.rm=TRUE),
    
    mean_neglog10FDR =
      mean(mean_neglog10FDR, na.rm=TRUE),
    
    n_sig_markers =
      mean(n_sig_markers, na.rm=TRUE),
    
    .groups="drop"
  )

# ----------------------
# Mean ± SD across datasets
# ----------------------
summary.stats <- summary.dataset %>%
  group_by(method) %>%
  summarise(
    
    logFC.mean = mean(median_log2FC),
    logFC.sd   = sd(median_log2FC),
    
    FDR.mean   = mean(mean_neglog10FDR),
    FDR.sd     = sd(mean_neglog10FDR),
    
    Mark.mean  = mean(n_sig_markers),
    Mark.sd    = sd(n_sig_markers)
    
  )

# ----------------------
# Colors
# ----------------------
method.cols <- c(
  Segmap="#FF0000",
  Scanpy="#008B00",
  Baysor="#8B2323",
  Seurat="#0000FF"
)

# ----------------------
# Panel A
# ----------------------
p1 <- ggplot() +
  geom_col(data=summary.stats,
           aes(method,FDR.mean, fill=method),
           width=.65, colour="black") +
  geom_errorbar(data=summary.stats,
                aes(method, ymin=FDR.mean-FDR.sd, ymax=FDR.mean+FDR.sd),
                width=.15) +
  
  geom_line(data=summary.dataset,
            aes(method, mean_neglog10FDR, group=Dataset),colour="grey70") +
  geom_point(data=summary.dataset,
             aes(method, mean_neglog10FDR), size=4.5, shape=21, fill="white") +
  scale_fill_manual(values=method.cols)+
  labs(title="Mean -log10(FDR)", x="Method",  y="-log10(FDR)")+
  theme_classic(base_size=15)+
  theme(legend.position="none",plot.title=element_text(face="bold"))
p1

# ----------------------
# Panel B
# ----------------------
p2 <- ggplot() +
  geom_col( data=summary.stats,
            aes(method, Mark.mean, fill=method),
            width=.65, colour="black") +
  geom_errorbar( data=summary.stats,
                 aes(method,
                     ymin=Mark.mean-Mark.sd,
                     ymax=Mark.mean+Mark.sd), width=.15) +
  geom_line(data=summary.dataset,
            aes(method, n_sig_markers, group=Dataset),
            colour="grey70") +
  geom_point( data=summary.dataset,
              aes(method, n_sig_markers),size=4.5,shape=21,fill="white") +
  scale_fill_manual(values=method.cols)+
  labs( title="Significant Marker Genes",
        x="Method", y="Number of markers") +
  theme_classic(base_size=15)+
  theme( legend.position="none", plot.title=element_text(face="bold"))
p2

# ----------------------
# Combine panels
# ----------------------
( p1 | p2) +
  plot_annotation(tag_levels="A")

# ==============================================================================
#   Figure 3c: Moran_per_gene
# ==============================================================================

data_Moran_per_gene_1 <- read.csv("segmap_data_Figure_3/data_Moran_per_gene_1.csv", header = TRUE, stringsAsFactors = FALSE)
data_Moran_per_gene_2 <- read.csv("segmap_data_Figure_3/data_Moran_per_gene_2.csv", header = TRUE, stringsAsFactors = FALSE)
data_Moran_per_gene_3 <- read.csv("segmap_data_Figure_3/data_Moran_per_gene_3.csv", header = TRUE, stringsAsFactors = FALSE)
data_Moran_per_gene_4 <- read.csv("segmap_data_Figure_3/data_Moran_per_gene_4.csv", header = TRUE, stringsAsFactors = FALSE)

# Make sure methods appear in desired order
data_Moran_per_gene_1$method <- factor(
  data_Moran_per_gene_1$method,
  levels = c("Segmap", "Scanpy", "Seurat", "Baysor")
)
data_Moran_per_gene_2$method <- factor(
  data_Moran_per_gene_2$method,
  levels = c("Segmap", "Scanpy", "Seurat", "Baysor")
)
data_Moran_per_gene_3$method <- factor(
  data_Moran_per_gene_3$method,
  levels = c("Segmap", "Scanpy", "Seurat", "Baysor")
)
data_Moran_per_gene_4$method <- factor(
  data_Moran_per_gene_4$method,
  levels = c("Segmap", "Scanpy", "Seurat", "Baysor")
)

data_Moran_per_gene_1$method[is.na(data_Moran_per_gene_1$method)] <- "Baysor"
data_Moran_per_gene_2$method[is.na(data_Moran_per_gene_2$method)] <- "Baysor"
data_Moran_per_gene_3$method[is.na(data_Moran_per_gene_3$method)] <- "Baysor"
data_Moran_per_gene_4$method[is.na(data_Moran_per_gene_4$method)] <- "Baysor"

ggplot(data_Moran_per_gene_1, aes(method, moran, fill = method)) +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = c(
    Segmap="#FF0000",
    Scanpy="#008B00",
    Baysor="#8B2323",
    Seurat="#0000FF"
  )) +
  labs(
    x = "Method",
    y = "Moran's I",
    title = "Dataset 1"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggplot(data_Moran_per_gene_2, aes(method, moran, fill = method)) +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = c(
    Segmap="#FF0000",
    Scanpy="#008B00",
    Baysor="#8B2323",
    Seurat="#0000FF"
  )) +
  labs(
    x = "Method",
    y = "Moran's I",
    title = "Dataset 2"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggplot(data_Moran_per_gene_3, aes(method, moran, fill = method)) +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = c(
    Segmap="#FF0000",
    Scanpy="#008B00",
    Baysor="#8B2323",
    Seurat="#0000FF"
  )) +
  labs(
    x = "Method",
    y = "Moran's I",
    title = "Dataset 3"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

ggplot(data_Moran_per_gene_4, aes(method, moran, fill = method)) +
  geom_boxplot(
    width = 0.55,
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.18,
    size = 0.7
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = c(
    Segmap="#FF0000",
    Scanpy="#008B00",
    Baysor="#8B2323",
    Seurat="#0000FF"
  )) +
  labs(
    x = "Method",
    y = "Moran's I",
    title = "Dataset 4"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

# ==============================================================================
#   Figure 3d: Spatial confidence-interval 
# ==============================================================================

# ----------------------
# Combine datasets
# ----------------------
data_spatial_CIs_1 <- read.csv("segmap_data_Figure_3/segmap_Benchmark_spatial_CIs_1.csv", header = TRUE, stringsAsFactors = FALSE)
data_spatial_CIs_2 <- read.csv("segmap_data_Figure_3/segmap_Benchmark_spatial_CIs_2.csv", header = TRUE, stringsAsFactors = FALSE)
data_spatial_CIs_3 <- read.csv("segmap_data_Figure_3/segmap_Benchmark_spatial_CIs_3.csv", header = TRUE, stringsAsFactors = FALSE)
data_spatial_CIs_4 <- read.csv("segmap_data_Figure_3/segmap_Benchmark_spatial_CIs_4.csv", header = TRUE, stringsAsFactors = FALSE)

data_spatial_CIs_1$Dataset <- "Dataset 1"
data_spatial_CIs_2$Dataset <- "Dataset 2"
data_spatial_CIs_3$Dataset <- "Dataset 3"
data_spatial_CIs_4$Dataset <- "Dataset 4"

all_data <- bind_rows(
  data_spatial_CIs_1,
  data_spatial_CIs_2,
  data_spatial_CIs_3,
  data_spatial_CIs_4
)

# ----------------------
# Keep only point estimates
# ----------------------
heatmap_data <- all_data %>%
  select(Dataset, Method, metric, point)

# ----------------------
# Rename metrics
# ----------------------
heatmap_data <- heatmap_data %>%
  mutate(
    metric = recode(
      metric,
      median_moran  = "Median Moran's I",
      pct_sig_moran = "% Sig. Moran",
      sp_magnitude  = "Spatial Magnitude"
    )
  )

# ----------------------
## Order factors
# ----------------------
heatmap_data$Dataset <- factor(
  heatmap_data$Dataset,
  levels = c(
    "Dataset 1",
    "Dataset 2",
    "Dataset 3",
    "Dataset 4"
  )
)

heatmap_data$Method <- factor(
  heatmap_data$Method,
  levels = c(
    "Scanpy",
    "Segmap",
    "Seurat",
    "Baysor"
  )
)

heatmap_data$metric <- factor(
  heatmap_data$metric,
  levels = c(
    "Median Moran's I",
    "% Sig. Moran",
    "Spatial Magnitude"
  )
)

# ----------------------
## Normalize within each metric and dataset
# ----------------------
heatmap_data <- heatmap_data %>%
  group_by(Dataset, metric) %>%
  mutate(
    score = (point - min(point)) /
      (max(point) - min(point))
  ) %>%
  ungroup()
heatmap_data <- na.omit(heatmap_data)

# ----------------------
# Plot
# ----------------------
p <-ggplot(heatmap_data, aes(metric, Method, fill = score)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_wrap( ~Dataset, nrow = 1) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0,1),
    name = "Normalized\nscore") +
  
  theme_classic(base_size = 16) +
  labs(
    x = "Spatial Metrics",
    y = "Segmentation Methods"
  ) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 18
    ),
    
    panel.border = element_blank(),
    panel.spacing = unit(0.4, "cm"),
    axis.title.x = element_text(
      # face = "bold",
      size = 16,
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      #face = "bold",
      size = 16,
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1,
      #face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      #face = "bold",
      size = 13
    ),
    legend.position = "right"
  )
p

# ==============================================================================
#                                Figue 3e: Bootstrap
# ==============================================================================

# ----------------------
# Combine datasets
# ----------------------
data_bootstrap_1 <- read.csv("segmap_data_Figure_3/data_bootstrap_1.csv", header = TRUE)
data_bootstrap_2 <- read.csv("segmap_data_Figure_3/data_bootstrap_2.csv", header = TRUE)
data_bootstrap_3 <- read.csv("segmap_data_Figure_3/data_bootstrap_3.csv", header = TRUE)
data_bootstrap_4 <- read.csv("segmap_data_Figure_3/data_bootstrap_4.csv", header = TRUE)

data_bootstrap_1$Dataset <- "Dataset 1"
data_bootstrap_2$Dataset <- "Dataset 2"
data_bootstrap_3$Dataset <- "Dataset 3"
data_bootstrap_4$Dataset <- "Dataset 4"

benchmark <- bind_rows(
  data_bootstrap_1,
  data_bootstrap_2,
  data_bootstrap_3,
  data_bootstrap_4
)

benchmark$Method <- factor(
  benchmark$Method,
  levels = c("Segmap","Scanpy","Seurat","Baysor")
)

# ----------------------
# Compute mean and SD
# ----------------------
summary.df <- benchmark %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(mean_rank),
    SD   = sd(mean_rank),
    .groups = "drop"
  )

summary.df

# ----------------------
# Define colors
# ----------------------
method.cols <- c(
  Segmap="#FF0000",
  Scanpy="#008B00",
  Baysor="#8B2323",
  Seurat="#0000FF"
)

summary.df <- benchmark %>%
  group_by(Method) %>%
  summarise(
    Mean = mean(mean_rank),
    SD = sd(mean_rank),
    .groups = "drop"
  )

ggplot() +
  geom_col( data = summary.df,
            aes(Method, Mean, fill = Method),
            width = 0.60,
            alpha = 0.75,
            color = "black") +
  
  geom_errorbar(data = summary.df,
                aes(Method, ymin = Mean - SD, ymax = Mean + SD), width = 0.15, linewidth = 0.8) +
  geom_line(data = benchmark,
            aes(Method, mean_rank, group = Dataset),
            colour = "grey70", linewidth = 0.6) +
  geom_point(data = benchmark,
             aes(Method, mean_rank),
             shape = 21,size = 4.5, stroke = 0.7, fill = "white", colour = "black") +
  scale_fill_manual(values = method.cols) +
  labs(x = NULL,
       y = "Mean Bootstrap Rank (Lower is Better)",
       title = "Benchmark Performance Across Four Datasets"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

# ==============================================================================
#               Figure 3f: weight sensitivity
# ==============================================================================
data_weight_sensitivity_1 <- read.csv("segmap_data_Figure_3/ddata_weight_sensitivity_1.csv", header = TRUE, stringsAsFactors = FALSE)
data_weight_sensitivity_2 <- read.csv("segmap_data_Figure_3/ddata_weight_sensitivity_2.csv", header = TRUE, stringsAsFactors = FALSE)
data_weight_sensitivity_3 <- read.csv("segmap_data_Figure_3/ddata_weight_sensitivity_3.csv", header = TRUE, stringsAsFactors = FALSE)
data_weight_sensitivity_4 <- read.csv("segmap_data_Figure_3/ddata_weight_sensitivity_4.csv", header = TRUE, stringsAsFactors = FALSE)

# ----------------------
# Combine datasets
# ----------------------
data_weight_sensitivity_1$Dataset <- "Dataset 1"
data_weight_sensitivity_2$Dataset <- "Dataset 2"
data_weight_sensitivity_3$Dataset <- "Dataset 3"
data_weight_sensitivity_4$Dataset <- "Dataset 4"

weight_data <-
  bind_rows(
    data_weight_sensitivity_1,
    data_weight_sensitivity_2,
    data_weight_sensitivity_3,
    data_weight_sensitivity_4
  )

# ----------------------
# Rename schemes
# ----------------------

weight_data <- weight_data %>%
  mutate(
    scheme = recode(
      scheme,
      as_specified = "As specified",
      equal        = "Equal",
      de_heavy     = "DE-heavy",
      sp_heavy     = "SP-heavy"
    )
  )

# ----------------------
# Order factors
# ----------------------
weight_data$Method <- factor(
  weight_data$Method,
  levels=c(
    "Scanpy",
    "Segmap",
    "Seurat"
  )
)

weight_data$scheme <- factor(
  weight_data$scheme,
  levels=c(
    "As specified",
    "Equal",
    "DE-heavy",
    "SP-heavy"
  )
)

weight_data$Dataset <- factor(
  weight_data$Dataset,
  levels=c(
    "Dataset 1",
    "Dataset 2",
    "Dataset 3",
    "Dataset 4"
  )
)

# ----------------------
# Heatmap
# ----------------------
ggplot(weight_data, aes(scheme, Method, fill=Overall))+
  geom_tile(colour="white",linewidth=0.6)+
  facet_wrap(~Dataset, nrow=1)+
  scale_fill_viridis_c(option="viridis", name="Overall\nScore")+
  labs(x="Weighting Scheme", y="Method") +
  theme_classic(base_size=16)+
  theme(strip.background=element_blank(),
        strip.text=element_text(face="bold", size=18),
        panel.border=element_blank(),
        panel.spacing = unit(0.4,"cm"),
        axis.text.x=element_text(
          angle=30,
          hjust=1,
          #face="bold",
          size=12),
        axis.text.y=element_text(
          #face="bold",
          size=13),
        
        axis.title=
          element_text(#face="bold",
            size=15),
        legend.position="right")

# ==============================================================================
#               Figure 3g: Per-method spatial coherence
# ==============================================================================

paired_bootstrap_compare <- function(
    moran_pergene,
    stat       = c("median_moran", "mean_moran", "pct_sig_moran", "sp_magnitude"),
    methods    = NULL,
    B          = 1000,
    conf       = 0.95,
    sp_weights = c(median_moran = 0.60, pct_sig_moran = 0.40),
    p_adjust   = "BH",
    seed       = 1
) {
  stat <- match.arg(stat)
  set.seed(seed)
  
  # ---- 1. Build the PAIRED gene x method matrices (genes present for ALL) ----
  if (is.null(methods)) methods <- sort(unique(moran_pergene$method))
  moran_pergene <- dplyr::filter(moran_pergene, method %in% methods)
  
  moran_wide <- moran_pergene %>%
    dplyr::select(gene, method, moran) %>%
    tidyr::pivot_wider(names_from = method, values_from = moran)
  pval_wide <- moran_pergene %>%
    dplyr::select(gene, method, pval) %>%
    tidyr::pivot_wider(names_from = method, values_from = pval)
  
  # keep only genes evaluated for EVERY method (valid pairing) ---------------
  keep <- stats::complete.cases(moran_wide[, methods, drop = FALSE])
  n_all <- nrow(moran_wide); n_keep <- sum(keep)
  if (n_keep < n_all)
    message(sprintf(paste0("paired_bootstrap_compare: %d of %d genes are present ",
                           "for all methods and used for the paired comparison (%d dropped)."),
                    n_keep, n_all, n_all - n_keep))
  if (n_keep < 10) stop("Fewer than 10 genes shared across all methods; paired ",
                        "comparison is not meaningful.")
  
  Mmoran <- as.matrix(moran_wide[keep, methods, drop = FALSE])   # genes x methods
  Mpval  <- as.matrix(pval_wide[keep,  methods, drop = FALSE])
  G <- nrow(Mmoran); K <- length(methods)
  
  # ---- 2. Per-method statistic from a (moran, pval) column ------------------
  sp_w <- sp_weights / sum(sp_weights)
  stat_fun <- switch(stat,
                     median_moran  = function(mo, pv) stats::median(mo, na.rm = TRUE),
                     mean_moran    = function(mo, pv) mean(mo, na.rm = TRUE),
                     pct_sig_moran = function(mo, pv) 100 * mean(p.adjust(pv, "BH") < 0.05, na.rm = TRUE),
                     sp_magnitude  = function(mo, pv) {
                       med  <- stats::median(mo, na.rm = TRUE)
                       psig <- mean(p.adjust(pv, "BH") < 0.05, na.rm = TRUE)   # fraction in [0,1]
                       as.numeric(sp_w["median_moran"] * med + sp_w["pct_sig_moran"] * psig)
                     })
  
  obs <- vapply(seq_len(K), function(k) stat_fun(Mmoran[, k], Mpval[, k]), numeric(1))
  names(obs) <- methods
  
  # ---- 3. Paired bootstrap: ONE gene resample shared across all methods -----
  boot <- matrix(NA_real_, B, K, dimnames = list(NULL, methods))
  for (b in seq_len(B)) {
    idx <- sample.int(G, G, replace = TRUE)          # shared indices = pairing
    mo_b <- Mmoran[idx, , drop = FALSE]
    pv_b <- Mpval[idx, , drop = FALSE]
    boot[b, ] <- vapply(seq_len(K), function(k) stat_fun(mo_b[, k], pv_b[, k]), numeric(1))
  }
  
  a <- (1 - conf) / 2
  per_method <- tibble(
    method  = methods,
    stat    = stat,
    observed = as.numeric(obs),
    ci_low  = apply(boot, 2, quantile, a,      na.rm = TRUE),
    ci_high = apply(boot, 2, quantile, 1 - a,  na.rm = TRUE)
  )
  
  # ---- 4. Pairwise DIFFERENCES from the paired bootstrap --------------------
  pairs <- utils::combn(methods, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    i <- p[1]; j <- p[2]
    d_obs  <- obs[[i]] - obs[[j]]
    d_boot <- boot[, i] - boot[, j]                  # paired difference dist.
    ci <- quantile(d_boot, c(a, 1 - a), na.rm = TRUE)
    # two-sided bootstrap p: mass on the side opposite the observed sign, x2,
    # floored at 1/B (resolution limit of B resamples).
    p_two <- 2 * min(mean(d_boot <= 0, na.rm = TRUE),
                     mean(d_boot >= 0, na.rm = TRUE))
    p_two <- min(max(p_two, 1 / B), 1)
    tibble(method_a = i, method_b = j,
           diff = d_obs, ci_low = ci[[1]], ci_high = ci[[2]],
           p_value = p_two)
  })
  pairwise <- dplyr::bind_rows(rows)
  pairwise$p_adj <- p.adjust(pairwise$p_value, method = p_adjust)
  pairwise <- pairwise %>%
    mutate(significant = (ci_low > 0) | (ci_high < 0),
           verdict = dplyr::case_when(
             significant &  diff > 0 ~ paste0(method_a, " > ", method_b),
             significant &  diff < 0 ~ paste0(method_b, " > ", method_a),
             TRUE                    ~ "tie (CI spans 0)")) %>%
    arrange(p_adj)
  
  list(pairwise = pairwise, per_method = per_method, boot = boot,
       n_genes = G, stat = stat)
}

# ============================
# CONFIG  -- the ONLY place you should need to edit
# ============================
CFG <- list(
  top_n         = 20,       # markers per cluster for the DE-strength summary
  moran_k       = 6,        # kNN for spatial weights
  n_boot        = 1000,     # bootstrap replicates
  fdr_floor     = 1e-300,   # floor so -log10(FDR) is finite
  
  # ---- DE-strength benchmark inputs (your "..._1" tables) ----
  de = list(
    seurat_rdata = "XeniumData_seurat/new_xenium_marker4.rdata",   # loads new_xenium_marker1
    seurat_obj   = "new_xenium_marker4",
    baysor_rdata = "XeniumData_baysor/baysor_marker4.rdata",
    baysor_obj = "baysor_marker4",
    scanpy_csv   = "XeniumData_scanpy/scanpy_data_4/marker_genes.csv",
    segmap_csv   = "XeniumData_segmap/cellpose_pipeline_outputs_4/data/cluster_markers_all_cells_final.csv"
  ),
  
  # ---- Spatial benchmark inputs (your "..._4" tables) ----
  sp = list(
    seurat_data_rdata   = "XeniumData_seurat/new_xenium_data4.rdata",   
    seurat_data_obj     = "new_xenium_data4",
    seurat_marker_rdata = "XeniumData_seurat/new_xenium_marker4.rdata", 
    seurat_marker_obj   = "new_xenium_marker4",
    
    baysor_data_rdata = "XeniumData_baysor/baysor_data_4.rdata",
    baysor_data_obj = "baysor_data_4",
    baysor_marker_rdata = "XeniumData_baysor/baysor_marker4.rdata", 
    baysor_marker_obj   = "baysor_marker4",
    
    scanpy_markers = "XeniumData_scanpy/scanpy_data_4/marker_genes.csv",
    scanpy_coords  = "XeniumData_scanpy/scanpy_data_4/spatial_coords.csv",
    scanpy_cellgene= "XeniumData_scanpy/scanpy_data_4/cell_by_gene.csv",
    segmap_markers = "XeniumData_segmap/cellpose_pipeline_outputs_4/data/cluster_markers_all_cells_final.csv",
    segmap_coords  = "XeniumData_segmap/cellpose_pipeline_outputs_4/data/spatial_clusters_final.csv",
    segmap_cellgene= "XeniumData_segmap/cellpose_pipeline_outputs_4/data/cell_by_gene_for_umap_final.csv"
  ),
  
  # ---- Column-name maps (edit if your headers differ) ----
  # pct1/pct2 = detection fraction inside / outside the cluster. Leave the
  # candidate vectors as-is; the first name that exists is used, else specificity
  # is treated as unavailable for that method.
  seurat_cols = list(cluster="cluster", gene="gene", log2FC="avg_log2FC",
                     fdr="p_val_adj", pct1=c("pct.1"), pct2=c("pct.2")),
  
  baysor_cols = list(cluster="cluster", gene="gene", log2FC="avg_log2FC",
                     fdr="p_val_adj", pct1=c("pct.1"), pct2=c("pct.2")),
  
  scanpy_cols = list(cluster="group", gene="names", log2FC="logfoldchanges",
                     fdr="pvals_adj",
                     pct1=c("pct_nz_group","pct.1","pct1","pct_in"),
                     pct2=c("pct_nz_reference","pct.2","pct2","pct_out")),
  
  segmap_cols = list(cluster="group", gene="names", log2FC="logfoldchanges",
                     fdr="pvals_adj",
                     pct1=c("pct_nz_group","pct.1","pct1","detection_in","pct_in"),
                     pct2=c("pct_nz_reference","pct.2","pct2","detection_out","pct_out")),
  
  # Coord column names per method
  seurat_coord_cols = c("x","y"),
  baysor_coord_cols = c("x","y"),
  scanpy_coord_cols = c("x","y"),
  segmap_coord_cols = c("x_centroid_px","y_centroid_px"),
  
  # Seurat assay/layer to read counts from for the SHARED normalisation.
  # We re-normalise all three the same way (CP10k + log1p) for Moran's I so the
  # comparison isn't confounded by SCT vs raw. If RNA counts are unavailable,
  # set seurat_moran_use_sct=TRUE to fall back to SCT 'data' (documented caveat).
  seurat_moran_assay = "RNA", seurat_moran_layer = "counts",
  seurat_moran_use_sct = FALSE,
  
  # ---- Metric weights (higher = better for every metric) ----
  # Within-benchmark weights, then across-benchmark weights. All are made
  # explicit and stress-tested in the sensitivity analysis below.
  w_de = c(median_log2FC = 0.40, neglog10FDR = 0.30,
           specificity   = 0.20, marker_yield = 0.10),
  w_sp = c(median_moran  = 0.60, pct_sig_moran = 0.40),
  w_overall = c(de = 0.5, sp = 0.5)   # DE strength vs spatial coherence
)

# ============================
# Small helpers
# ============================
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Pick the first candidate column that exists in df; else NA_character_
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit)) hit[1] else NA_character_
}

# CP10k + log1p on a genes x cells (sparse) count matrix -> genes x cells
lognorm_counts <- function(counts) {
  counts <- as(counts, "CsparseMatrix")
  lib <- Matrix::colSums(counts); lib[lib == 0] <- 1
  cp10k <- counts %*% Matrix::Diagonal(x = 1e4 / lib)
  # log1p keeps zeros zero -> stays sparse
  cp10k@x <- log1p(cp10k@x)
  cp10k
}

# rank-score: best -> 1, worst -> 0 (NA preserved). higher_better as named.
rank_score <- function(x, higher_better = TRUE) {
  x <- as.numeric(x)
  ok <- !is.na(x); n <- sum(ok)
  out <- rep(NA_real_, length(x))
  if (n <= 1) { out[ok] <- if (n == 1) 1 else NA_real_; return(out) }
  r <- rank(if (higher_better) -x[ok] else x[ok], ties.method = "average")
  out[ok] <- (n - r) / (n - 1)
  out
}

# integer competition rank (1 = best), higher_better
rank_int <- function(x, higher_better = TRUE) {
  rank(if (higher_better) -x else x, ties.method = "min", na.last = "keep")
}

# ============================
# PART 1 -- Load & standardise marker tables into one tidy schema
#   -> columns: method, cluster, gene, log2FC, fdr, pct1, pct2  (pct* may be NA)
# ============================
load_rdata_obj <- function(path, objname) {
  e <- new.env(); load(path, envir = e)
  if (!is.null(objname) && objname %in% ls(e)) return(get(objname, envir = e))
  get(ls(e)[1], envir = e)  # fall back to the single object in the file
}

standardise_markers <- function(df, method, cols) {
  df <- as.data.frame(df)
  p1 <- pick_col(df, cols$pct1); p2 <- pick_col(df, cols$pct2)
  spec_ok <- !is.na(p1) && !is.na(p2)
  tibble(
    method  = method,
    cluster = as.character(df[[cols$cluster]]),
    gene    = as.character(df[[cols$gene]]),
    log2FC  = suppressWarnings(as.numeric(df[[cols$log2FC]])),
    fdr     = suppressWarnings(as.numeric(df[[cols$fdr]])),
    pct1    = if (spec_ok) suppressWarnings(as.numeric(df[[p1]])) else NA_real_,
    pct2    = if (spec_ok) suppressWarnings(as.numeric(df[[p2]])) else NA_real_
  )
}

message("Loading marker tables (DE-strength benchmark) ...")
mk_seurat_raw <- load_rdata_obj(CFG$de$seurat_rdata, CFG$de$seurat_obj)
mk_baysor_raw <- load_rdata_obj(CFG$de$baysor_rdata, CFG$de$baysor_obj)
mk_scanpy_raw <- readr::read_csv(CFG$de$scanpy_csv, show_col_types = FALSE)
mk_segmap_raw <- readr::read_csv(CFG$de$segmap_csv, show_col_types = FALSE)

markers <- bind_rows(
  standardise_markers(mk_seurat_raw, "Seurat", CFG$seurat_cols),
  standardise_markers(mk_baysor_raw, "Baysor", CFG$baysor_cols),
  standardise_markers(mk_scanpy_raw, "Scanpy", CFG$scanpy_cols),
  standardise_markers(mk_segmap_raw, "Segmap", CFG$segmap_cols)
)


# Which methods actually have specificity? (drives the "drop metric for all" rule)
spec_available <- markers %>% group_by(method) %>%
  summarise(has_spec = any(!is.na(pct1) & !is.na(pct2)), .groups = "drop")
USE_SPECIFICITY <- all(spec_available$has_spec)
if (!USE_SPECIFICITY) {
  message("NOTE: specificity (pct.1-pct.2) is not available for all methods (",
          paste(spec_available$method[!spec_available$has_spec], collapse = ", "),
          "). Dropping specificity from the composite for ALL methods to keep the ",
          "comparison neutral. To include it, add detection-fraction columns to ",
          "those marker exports (Scanpy: rank_genes_groups(..., pts=TRUE)).")
}

# ============================
# PART 2 -- DE-strength: per-cluster summary + method-level metrics
# ============================

# significant POSITIVE markers only
sig <- markers %>% filter(is.finite(fdr), fdr < 0.05, is.finite(log2FC), log2FC > 0)

# marker YIELD is counted BEFORE the top-N cap (the cap hides the real difference)
yield_percluster <- sig %>% count(method, cluster, name = "n_sig_markers")

# top-N per cluster for logFC / FDR / specificity summaries
topN <- sig %>% group_by(method, cluster) %>%
  slice_max(order_by = log2FC, n = CFG$top_n, with_ties = FALSE) %>% ungroup()

de_percluster <- topN %>%
  mutate(neglog10FDR = -log10(pmax(fdr, CFG$fdr_floor)),
         specificity = pct1 - pct2) %>%
  group_by(method, cluster) %>%
  summarise(median_log2FC = median(log2FC),
            mean_neglog10FDR = mean(neglog10FDR),
            specificity = if (USE_SPECIFICITY) mean(specificity) else NA_real_,
            n_top = n(), .groups = "drop") %>%
  left_join(yield_percluster, by = c("method", "cluster")) %>%
  mutate(n_sig_markers = tidyr::replace_na(n_sig_markers, 0L))

readr::write_csv(de_percluster, "segmap_Benchmark_DE_per_cluster.csv")

# ============================
# PART 3 -- Spatial coherence (Moran's I) on a COMMON gene set + COMMON norm
# ============================

# Build a listw from a coordinate matrix
build_listw <- function(coords, k) {
  nb <- spdep::knn2nb(spdep::knearneigh(as.matrix(coords), k = k))
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}

# Moran's I for one gene vector given a listw (returns c(I, p) or NA)
moran_one <- function(x, lw) {
  if (anyNA(x) || stats::var(x) == 0 || sum(x > 0) < 10) return(c(NA_real_, NA_real_))
  res <- tryCatch(spdep::moran.test(x, lw, zero.policy = TRUE),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(res)) return(c(NA_real_, NA_real_))
  c(unname(res$estimate["Moran I statistic"]), res$p.value)
}

# Compute Moran's I for a set of genes on a (genes x cells) normalised matrix
moran_over_genes <- function(expr_norm, lw, genes, method) {
  genes <- genes[genes %in% rownames(expr_norm)]
  res <- lapply(genes, function(g) {
    v <- moran_one(as.numeric(expr_norm[g, ]), lw)
    tibble(method = method, gene = g, moran = v[1], pval = v[2])
  })
  bind_rows(res) %>% filter(!is.na(moran))
}

message("Loading spatial inputs (Moran's I benchmark) ...")

## ----- Seurat: object -> coords + normalised expression -----
seurat_obj <- load_rdata_obj(CFG$sp$seurat_data_rdata, CFG$sp$seurat_data_obj)
seurat_coords <- as.matrix(GetTissueCoordinates(seurat_obj)[, CFG$seurat_coord_cols])
lw_seurat <- build_listw(seurat_coords, CFG$moran_k)

if (!CFG$seurat_moran_use_sct) {
  expr_seurat <- tryCatch(
    lognorm_counts(GetAssayData(seurat_obj, assay = CFG$seurat_moran_assay,
                                layer = CFG$seurat_moran_layer)),
    error = function(e) {
      message("  Seurat RNA counts not found; falling back to SCT 'data' ",
              "(normalisation not identical to Scanpy/Segmap).")
      GetAssayData(seurat_obj, assay = "SCT", layer = "data")
    })
} else {
  expr_seurat <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
}

## ----- Baysor: object -> coords + normalised expression -----
baysor_obj <- load_rdata_obj(CFG$sp$baysor_data_rdata, CFG$sp$baysor_data_obj)
baysor_coords <- as.matrix(baysor_obj@meta.data[, c("x", "y")])
lw_baysor <- build_listw(baysor_coords, CFG$moran_k)

# BAYSOR MORAN CONFIGURATION
CFG$baysor_moran_use_sct <- "SCT" %in% Assays(baysor_obj)

if (CFG$baysor_moran_use_sct) {
  CFG$baysor_moran_assay <- "SCT"
  CFG$baysor_moran_layer <- "data"
} else {
  CFG$baysor_moran_assay <- "RNA"
  CFG$baysor_moran_layer <- "counts"
}

# EXTRACT EXPRESSION MATRIX
if (CFG$baysor_moran_use_sct) {
  
  message("Using SCT assay.")
  
  expr_baysor <- GetAssayData(
    baysor_obj,
    assay = "SCT",
    layer = "data"
  )
  
} else {
  
  message("Using RNA counts.")
  
  expr_baysor <- tryCatch(
    
    lognorm_counts(
      GetAssayData(
        baysor_obj,
        assay = "RNA",
        layer = "counts"
      )
    ),
    
    error = function(e) {
      
      message("RNA counts layer not found. Using RNA data layer instead.")
      
      GetAssayData(
        baysor_obj,
        assay = "RNA",
        layer = "data"
      )
      
    }
    
  )
  
}

cat("Assay used :", CFG$baysor_moran_assay, "\n")
cat("Layer used :", CFG$baysor_moran_layer, "\n")
cat("Expression matrix dimensions:\n")
print(dim(expr_baysor))


## ----- Scanpy: coords + cell_by_gene -----
sc_coords <- readr::read_csv(CFG$sp$scanpy_coords, show_col_types = FALSE)
sc_cg     <- readr::read_csv(CFG$sp$scanpy_cellgene, show_col_types = FALSE)
lw_scanpy <- build_listw(as.matrix(sc_coords[, CFG$scanpy_coord_cols]), CFG$moran_k)
sc_counts <- Matrix::Matrix(t(as.matrix(sc_cg[, -1])), sparse = TRUE)
rownames(sc_counts) <- colnames(sc_cg)[-1]; colnames(sc_counts) <- sc_cg[[1]]
expr_scanpy <- lognorm_counts(sc_counts)

## ----- Segmap: coords + cell_by_gene -----
sg_coords <- readr::read_csv(CFG$sp$segmap_coords, show_col_types = FALSE)
sg_cg     <- readr::read_csv(CFG$sp$segmap_cellgene, show_col_types = FALSE)
lw_segmap <- build_listw(as.matrix(sg_coords[, CFG$segmap_coord_cols]), CFG$moran_k)
sg_counts <- Matrix::Matrix(t(as.matrix(sg_cg[, -1])), sparse = TRUE)
rownames(sg_counts) <- colnames(sg_cg)[-1]; colnames(sg_counts) <- sg_cg[[1]]
expr_segmap <- lognorm_counts(sg_counts)

## ----- COMMON gene set: genes present in all three matrices -----
common_genes <- Reduce(intersect, list(rownames(expr_seurat),
                                       rownames(expr_baysor),
                                       rownames(expr_scanpy),
                                       rownames(expr_segmap)))
message("Common genes evaluated for Moran's I (all methods, same list): ",
        length(common_genes))
stopifnot(length(common_genes) >= 10)

moran_pergene <- bind_rows(
  moran_over_genes(expr_seurat, lw_seurat, common_genes, "Seurat"),
  moran_over_genes(expr_baysor, lw_baysor, common_genes, "Baysor"),
  moran_over_genes(expr_scanpy, lw_scanpy, common_genes, "Scanpy"),
  moran_over_genes(expr_segmap, lw_segmap, common_genes, "Segmap")
)
res <- paired_bootstrap_compare(moran_pergene, stat = "median_moran", B = CFG$n_boot)


# ==============================================================================
#               Figure 3g: Per-method spatial coherence
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

# internal: human-readable label for the statistic stored in res$stat
.stat_label <- function(stat) {
  switch(stat,
         median_moran  = "median Moran's I",
         mean_moran    = "mean Moran's I",
         pct_sig_moran = "% genes with significant Moran's I",
         sp_magnitude  = "spatial magnitude",
         stat)
}

.default_method_colors <- c(Segmap = "#FF0000", Scanpy = "#008B00",
                            Seurat = "#0000FF", Baysor = "#8B2323")

# ---------------------------------------------------------------------------
# p1: per-method point estimate + 95% bootstrap CI (standalone)
# ---------------------------------------------------------------------------
plot_paired_per_method <- function(res,
                                   method_colors = .default_method_colors,
                                   digits = 4,
                                   title = "Per-method spatial coherence") {
  stopifnot("per_method" %in% names(res))
  stat_lab <- .stat_label(res$stat)
  
  pm <- res$per_method %>%
    mutate(method = factor(method, levels = method[order(observed)]))
  cols <- method_colors[as.character(levels(pm$method))]
  if (any(is.na(cols))) cols[is.na(cols)] <- "#666666"
  
  # --------------------------
  # User-defined plot parameters
  # --------------------------
  
  point_size      <- 5
  errorbar_size   <- 2
  text_size       <- 4
  title_size      <- 16
  axis_title_size <- 14
  axis_text_size  <- 12
  axis_line_size  <- 0.8
  
  x_label <- paste0(stat_lab, " (95% bootstrap CI)")
  y_label <- "Method"
  
  # --------------------------
  # Plot
  # --------------------------
  
  ggplot(pm, aes(x = observed, y = method, color = method)) +
    
    geom_errorbarh(
      aes(xmin = ci_low, xmax = ci_high),
      height = 0,
      linewidth = errorbar_size,
      na.rm = TRUE
    ) +
    
    geom_point(
      size = point_size,
      na.rm = TRUE
    ) +
    
    geom_text(
      aes(
        x = ci_high,
        label = formatC(
          observed,
          format = "f",
          digits = digits
        )
      ),
      hjust = -0.25,
      size = text_size,
      color = "grey25"
    ) +
    
    scale_color_manual(
      values = cols,
      guide = "none"
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.05, 0.22))
    ) +
    
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    
    theme_classic(base_size = axis_text_size) +
    
    theme(
      
      # Completely remove all grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # White backgrounds everywhere
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      legend.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      legend.key = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = axis_title_size,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = axis_title_size,
        face = "bold"
      ),
      
      # Axis text
      axis.text.x = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      
      # Axis lines
      axis.line = element_line(
        linewidth = axis_line_size,
        colour = "black"
      ),
      
      # Plot title
      plot.title = element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5
      )
    )
}

# ==============================================================================
#  Figure 3h: pairwise differences + 95% paired-bootstrap CI + zero line (standalone)
# ==============================================================================

plot_paired_differences <- function(res,
                                    digits = 2,
                                    title = "Direct paired comparisons  (CI excludes 0 = significant)") {
  stopifnot("pairwise" %in% names(res))
  stat_lab <- .stat_label(res$stat)
  
  pw <- res$pairwise %>%
    mutate(pair  = paste(method_a, "vs", method_b),
           pair  = factor(pair, levels = pair[order(diff)]),
           annot = paste0(verdict, "  (p=", formatC(p_adj, format = "g", digits = digits), ")"))
  
  # --------------------------
  # User-defined plot parameters
  # --------------------------
  
  point_size       <- 4
  errorbar_size    <- 2.5
  text_size        <- 4
  title_size       <- 18
  axis_title_size  <- 16
  axis_text_size   <- 14
  axis_line_size   <- 0.8
  vline_size       <- 2.5
  
  x_label <- paste0(
    "Difference in ", stat_lab,
    " (95% paired-bootstrap CI)"
  )
  
  y_label <- "Comparison"
  
  # --------------------------
  # Plot
  # --------------------------
  
  ggplot(pw, aes(x = diff, y = pair)) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#B0483B",
      linewidth = vline_size
    ) +
    
    geom_errorbarh(
      aes(xmin = ci_low, xmax = ci_high),
      height = 0,
      linewidth = errorbar_size,
      color = "#68228B",
      na.rm = TRUE
    ) +
    
    geom_point(
      shape = 18,
      size = point_size,
      color = "grey5",
      na.rm = TRUE
    ) +
    
    geom_text(
      aes(
        x = ci_high,
        label = annot
      ),
      hjust = -0.08,
      size = text_size,
      color = "grey20"
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0.08, 0.55))
    ) +
    
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    
    theme_classic(base_size = axis_text_size) +
    
    theme(
      
      # Completely remove grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # White backgrounds
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      plot.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      legend.background = element_rect(
        fill = "white",
        colour = NA
      ),
      
      legend.key = element_rect(
        fill = "white",
        colour = NA
      ),
      
      # Axis titles
      axis.title.x = element_text(
        size = axis_title_size,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = axis_title_size,
        face = "bold"
      ),
      
      # Axis labels
      axis.text.x = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      
      # Axis lines
      axis.line = element_line(
        linewidth = axis_line_size,
        colour = "black"
      ),
      
      # Title
      plot.title = element_text(
        size = title_size,
        face = "bold",
        hjust = 0.5
      )
    )
}

# ---------------------------------------------------------------------------
# combined helper (optional): p1 | p2 side by side
# ---------------------------------------------------------------------------
plot_paired_bootstrap <- function(res,
                                  method_colors = .default_method_colors,
                                  digits = 4) {
  p1 <- plot_paired_per_method(res, method_colors = method_colors, digits = digits)
  p2 <- plot_paired_differences(res)
  stat_lab <- .stat_label(res$stat)
  ttl <- sprintf("Direct paired bootstrap comparison of spatial coherence  (%s, %d common genes)",
                 stat_lab, res$n_genes)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, widths = c(1, 1.25)) +
      patchwork::plot_annotation(
        title = ttl,
        theme = theme(plot.title = element_text(face = "bold", size = 13)))
  } else {
    message("patchwork not installed; returning list(per_method = p1, pairwise = p2).")
    list(per_method = p1, pairwise = p2)
  }
}

plot_paired_per_method(res)
plot_paired_differences(res)


