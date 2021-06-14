#######################################################
# Script to generate figures for cytometry review paper
# Lukas Weber, June 2021
#######################################################

library(HDCytoData)
library(SummarizedExperiment)
library(FlowSOM)
library(Rtsne)
library(umap)
library(ggplot2)


# ---------
# Load data
# ---------

d_se <- Levine_32dim_SE()


# -------------
# Preprocessing
# -------------

# select cell type marker columns
d_clus <- assay(d_se)[, colData(d_se)$marker_class == "type"]
dim(d_clus)
d_clus[1:6, 1:6]

pop_ids <- rowData(d_se)$population_id
stopifnot(length(pop_ids) == nrow(d_clus))

# transform data
cofactor <- 5
d_clus <- asinh(d_clus / cofactor)
summary(d_clus)

# subsample for figures
n <- 10000
set.seed(123)
ix <- sample(seq_len(nrow(d_clus)), n)
d_sub <- d_clus[ix, ]

pop_ids_sub <- pop_ids[ix]
stopifnot(length(pop_ids_sub) == nrow(d_sub))

# remove any near-duplicate rows (required by Rtsne)
dups <- duplicated(d_sub)
table(dups)
d_sub <- d_sub[!dups, ]

pop_ids_sub <- pop_ids_sub[!dups]
stopifnot(length(pop_ids_sub) == nrow(d_sub))

ix_keep <- ix[!dups]


# -------------------
# Dimension reduction
# -------------------

n_dims <- 2


# run Rtsne with different random seeds

set.seed(101)
out_Rtsne_1 <- Rtsne(as.matrix(d_sub), dims = n_dims)
dims_Rtsne_1 <- out_Rtsne_1$Y
colnames(dims_Rtsne_1) <- c("tSNE_1", "tSNE_2")
head(dims_Rtsne_1)

stopifnot(nrow(dims_Rtsne_1) == length(pop_ids_sub))

set.seed(102)
out_Rtsne_2 <- Rtsne(as.matrix(d_sub), dims = n_dims)
dims_Rtsne_2 <- out_Rtsne_2$Y
colnames(dims_Rtsne_2) <- c("tSNE_1", "tSNE_2")
head(dims_Rtsne_2)

stopifnot(nrow(dims_Rtsne_2) == length(pop_ids_sub))

set.seed(103)
out_Rtsne_3 <- Rtsne(as.matrix(d_sub), dims = n_dims)
dims_Rtsne_3 <- out_Rtsne_3$Y
colnames(dims_Rtsne_3) <- c("tSNE_1", "tSNE_2")
head(dims_Rtsne_3)

stopifnot(nrow(dims_Rtsne_3) == length(pop_ids_sub))


# run UMAP with different random seeds

set.seed(101)
out_umap_1 <- umap(d_sub)
dims_umap_1 <- out_umap_1$layout
colnames(dims_umap_1) <- c("UMAP_1", "UMAP_2")
head(dims_umap_1)

stopifnot(nrow(dims_umap_1) == length(pop_ids_sub))

set.seed(102)
out_umap_2 <- umap(d_sub)
dims_umap_2 <- out_umap_2$layout
colnames(dims_umap_2) <- c("UMAP_1", "UMAP_2")
head(dims_umap_2)

stopifnot(nrow(dims_umap_2) == length(pop_ids_sub))

set.seed(103)
out_umap_3 <- umap(d_sub)
dims_umap_3 <- out_umap_3$layout
colnames(dims_umap_3) <- c("UMAP_1", "UMAP_2")
head(dims_umap_3)

stopifnot(nrow(dims_umap_3) == length(pop_ids_sub))


# ----------
# Clustering
# ----------

# clustering using FlowSOM

d_flowsom <- flowFrame(d_sub)

set.seed(123)

out_fs <- ReadInput(d_flowsom, transform = FALSE, scale = FALSE)
out_fs <- BuildSOM(out_fs)
out_fs <- BuildMST(out_fs)

labels_pre_meta <- out_fs$map$mapping[, 1]

# number of meta-clusters
k <- 12

seed <- 123
out_fs <- metaClustering_consensus(out_fs$map$codes, k = k, seed = seed)

labels_fs <- out_fs[labels_pre_meta]

# summary of clusters
table(labels_fs)


# --------------------------
# Generate plots: clustering
# --------------------------

library(RColorBrewer)
colors_clus <- brewer.pal(12, "Paired")

#library(colorspace)
#set.seed(7)
#colors_clus <- sample(qualitative_hcl(12))

#colors_clus <- unname(palette.colors(palette = "Polychrome 36"))
#colors_clus <- rainbow(12)
#colors_clus <- brewer.pal(12, "Set3")

# function for tSNE plots
plot_tSNE_clus <- function(dims_tSNE, labels, colors, 
                           size = 0.5, alpha = 0.5) {
  
  d_plot <- cbind(as.data.frame(dims_tSNE), cluster = as.factor(labels))
  
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
    geom_point(size = size, alpha = alpha) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(legend.key.size = unit(0.75, "lines"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))
}

# seed 1
plot_tSNE_clus(dims_Rtsne_1, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed1.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_tSNE_seed1.pdf", width = 4, height = 3)

# seed 2
plot_tSNE_clus(dims_Rtsne_2, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed2.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_tSNE_seed2.pdf", width = 4, height = 3)

# seed 3
plot_tSNE_clus(dims_Rtsne_3, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed3.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_tSNE_seed3.pdf", width = 4, height = 3)


# function for UMAP plots
plot_UMAP_clus <- function(dims_umap, labels, colors, 
                           size = 0.2, alpha = 0.5) {
  
  d_plot <- cbind(as.data.frame(dims_umap), cluster = as.factor(labels))
  
  ggplot(d_plot, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
    geom_point(size = size, alpha = alpha) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(legend.key.size = unit(0.75, "lines"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))
}

# seed 1
plot_UMAP_clus(dims_umap_1, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed1.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_UMAP_seed1.pdf", width = 4, height = 3)

# seed 2
plot_UMAP_clus(dims_umap_2, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed2.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_UMAP_seed2.pdf", width = 4, height = 3)

# seed 3
plot_UMAP_clus(dims_umap_3, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed3.png", width = 4, height = 3)
ggsave("../plots/FlowSOM_UMAP_seed3.pdf", width = 4, height = 3)


# ---------------
# Expression plot
# ---------------

# function to check expression of a given marker
plot_UMAP_expr <- function(dims_umap, d_sub, marker, 
                           size = 0.1, alpha = 0.5) {
  
  d_plot <- cbind(as.data.frame(dims_umap), expr = d_sub[, marker])
  
  ggplot(d_plot, aes(x = UMAP_1, y = UMAP_2, color = expr)) + 
    geom_point(size = size, alpha = alpha) + 
    scale_color_viridis_c() + 
    ggtitle(marker) + 
    theme_bw() + 
    theme(legend.key.size = unit(0.75, "lines"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

marker <- "CD14"
plot_UMAP_expr(dims_umap_1, d_sub, marker)
ggsave(paste0("../plots/", marker, "_UMAP.png"), width = 3.5, height = 3)

marker <- "CD16"
plot_UMAP_expr(dims_umap_1, d_sub, marker)
ggsave(paste0("../plots/", marker, "_UMAP.png"), width = 3.5, height = 3)


# ---------------------------
# Facetted plot: random seeds
# ---------------------------

# facetted plot for 2 seeds
d_umap <- rbind(
  cbind(as.data.frame(dims_umap_1), seed = "seed1", cluster = labels_fs), 
  cbind(as.data.frame(dims_umap_2), seed = "seed2", cluster = labels_fs))

d_umap$seed <- as.factor(d_umap$seed)
d_umap$cluster <- as.factor(d_umap$cluster)

ggplot(d_umap, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = colors_clus) + 
  facet_wrap(~seed, ncol = 1) + 
  theme_bw() + 
  theme(legend.key.size = unit(0.75, "lines"), 
        panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))

ggsave("../plots/UMAP_randomseeds_two.png", width = 3.75, height = 5.5)
ggsave("../plots/UMAP_randomseeds_two.pdf", width = 3.75, height = 5.5)


# facetted plot for 3 seeds
d_umap <- rbind(
  cbind(as.data.frame(dims_umap_1), seed = "seed1", cluster = labels_fs), 
  cbind(as.data.frame(dims_umap_2), seed = "seed2", cluster = labels_fs), 
  cbind(as.data.frame(dims_umap_3), seed = "seed3", cluster = labels_fs))

d_umap$seed <- as.factor(d_umap$seed)
d_umap$cluster <- as.factor(d_umap$cluster)

ggplot(d_umap, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = colors_clus) + 
  facet_wrap(~seed, ncol = 1) + 
  theme_bw() + 
  theme(legend.key.size = unit(0.75, "lines"), 
        panel.grid = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))

ggsave("../plots/UMAP_randomseeds_three.png", width = 3.5, height = 7.5)
ggsave("../plots/UMAP_randomseeds_three.pdf", width = 3.5, height = 7.5)


# ----------------------------
# Generate plots: ground truth
# ----------------------------

colors_truth <- c(rainbow(14), "gray85")

# function for tSNE plots
plot_tSNE_truth <- function(dims_tSNE, pop_ids, colors) {
  
  d_plot <- cbind(as.data.frame(dims_tSNE), population = as.factor(pop_ids))
  
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, color = population)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(legend.key.size = unit(0.75, "lines"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))
}

# seed 1
plot_tSNE_truth(dims_Rtsne_1, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed1.png", width = 6, height = 3.5)
ggsave("../plots/truth_tSNE_seed1.pdf", width = 6, height = 3.5)

# seed 2
plot_tSNE_truth(dims_Rtsne_2, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed2.png", width = 6, height = 3.5)
ggsave("../plots/truth_tSNE_seed2.pdf", width = 6, height = 3.5)

# seed 3
plot_tSNE_truth(dims_Rtsne_3, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed3.png", width = 6, height = 3.5)
ggsave("../plots/truth_tSNE_seed3.pdf", width = 6, height = 3.5)


# function for UMAP plots
plot_UMAP_truth <- function(dims_umap, pop_ids, colors) {
  
  d_plot <- cbind(as.data.frame(dims_umap), population = as.factor(pop_ids))
  
  ggplot(d_plot, aes(x = UMAP_1, y = UMAP_2, color = population)) + 
    geom_point(size = 0.2, alpha = 0.5) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(legend.key.size = unit(0.75, "lines"), 
          panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))
}

# seed 1
plot_UMAP_truth(dims_umap_1, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed1.png", width = 6, height = 3.5)
ggsave("../plots/truth_UMAP_seed1.pdf", width = 6, height = 3.5)

# seed 2
plot_UMAP_truth(dims_umap_2, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed2.png", width = 6, height = 3.5)
ggsave("../plots/truth_UMAP_seed2.pdf", width = 6, height = 3.5)

# seed 3
plot_UMAP_truth(dims_umap_3, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed3.png", width = 6, height = 3.5)
ggsave("../plots/truth_UMAP_seed3.pdf", width = 6, height = 3.5)


# --------------------------
# Cluster expression heatmap
# --------------------------

# use heatmap plotting function from CATALYST

library(CATALYST)
library(scales)

# format into SingleCellExperiment for CATALYST
sce <- SingleCellExperiment(
  rowData = colData(d_se)[colData(d_se)$marker_class == "type", ], 
  colData = rowData(d_se)[ix_keep, ], 
  assays = list(
    exprs = t(d_sub)
  )
)

colData(sce)$cluster_id <- labels_fs

# structure required by plotExprHeatmap
colData(sce)$sample_id <- "sample01"
metadata(sce)$cluster_codes <- data.frame(som100 = factor(1:k), meta12 = factor(1:k))

# order markers alphabetically
sce <- sce[order(rowData(sce)$marker_name), ]

pdf("../plots/heatmap_clusters.pdf", width = 6.5, height = 3.5)
plotExprHeatmap(sce, features = "type", by = "cluster_id", 
                k = "meta12", k_pal = colors_clus, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()


# alternative code

# library(dplyr)
# library(tidyr)
# 
# df <- cbind(as.data.frame(d_sub), cluster_id = labels_fs, cell_id = 1:nrow(d_sub)) %>% 
#   gather("marker", "expr", -cluster_id, -cell_id) %>% 
#   group_by(cluster_id, marker) %>% 
#   summarize(median = median(expr)) %>% 
#   spread(marker, median)
# 
# mat <- as.matrix(df[, 2:33])
# rownames(mat) <- df$cluster_id


# -------------------------
# Clustering within T cells
# -------------------------

# identify T cells from FlowSOM clustering by inspecting UMAP and heatmap

# T cells are clusters 1 and 2
table(labels_fs)
# number and proportion of cells identified as T cells
table(labels_fs %in% c(1, 2))
mean(labels_fs %in% c(1, 2))

d_sub_tcells <- d_sub[labels_fs %in% c(1, 2), ]
dim(d_sub_tcells)

# perform second round of clustering on T cells only
d_flowsom_tcells <- flowFrame(d_sub_tcells)
set.seed(123)
out_fs_tcells <- ReadInput(d_flowsom_tcells, transform = FALSE, scale = FALSE)
out_fs_tcells <- BuildSOM(out_fs_tcells)
out_fs_tcells <- BuildMST(out_fs_tcells)
labels_pre_meta_tcells <- out_fs_tcells$map$mapping[, 1]
# number of meta-clusters
k <- 12
seed <- 123
out_fs_tcells <- metaClustering_consensus(out_fs_tcells$map$codes, k = k, seed = seed)
labels_fs_tcells <- out_fs_tcells[labels_pre_meta_tcells]
# summary of clusters
table(labels_fs_tcells)


# calculate tSNE and UMAP dimension reduction: T cells

n_dims <- 2

set.seed(100)
out_Rtsne <- Rtsne(as.matrix(d_sub_tcells), dims = n_dims)
dims_Rtsne <- out_Rtsne$Y
colnames(dims_Rtsne) <- c("tSNE_1", "tSNE_2")
stopifnot(nrow(dims_Rtsne) == length(labels_fs_tcells))

set.seed(100)
out_umap_tcells <- umap(d_sub_tcells)
dims_umap_tcells <- out_umap_tcells$layout
colnames(dims_umap_tcells) <- c("UMAP_1", "UMAP_2")
stopifnot(nrow(dims_umap_tcells) == length(labels_fs_tcells))


# plot T cell subsets: tSNE and UMAP

colors_tcells <- qualitative_hcl(k)

plot_tSNE_clus(dims_Rtsne, labels_fs_tcells, colors_tcells)
ggsave("../plots/Tcells_FlowSOM_tSNE.png", width = 4, height = 3)
ggsave("../plots/Tcells_FlowSOM_tSNE.pdf", width = 4, height = 3)

plot_UMAP_clus(dims_umap_tcells, labels_fs_tcells, colors_tcells)
ggsave("../plots/Tcells_FlowSOM_UMAP.png", width = 4, height = 3)
ggsave("../plots/Tcells_FlowSOM_UMAP.pdf", width = 4, height = 3)


# heatmap: T cells

ix_tcells <- labels_fs %in% c(1, 2)
stopifnot(length(ix_tcells) == length(ix_keep))
stopifnot(nrow(d_sub_tcells) == nrow(rowData(d_se)[ix_keep, ][ix_tcells, ]))

# format into SingleCellExperiment for CATALYST
sce_tcells <- SingleCellExperiment(
  rowData = colData(d_se)[colData(d_se)$marker_class == "type", ], 
  colData = rowData(d_se)[ix_keep, ][ix_tcells, ], 
  assays = list(
    exprs = t(d_sub_tcells)
  )
)

colData(sce_tcells)$cluster_id <- labels_fs_tcells

# structure required by plotExprHeatmap
colData(sce_tcells)$sample_id <- "sample01"
metadata(sce_tcells)$cluster_codes <- data.frame(som100 = factor(1:k), meta12 = factor(1:k))

# order markers alphabetically
sce_tcells <- sce_tcells[order(rowData(sce_tcells)$marker_name), ]

pdf("../plots/heatmap_tcells.pdf", width = 6.5, height = 3.5)
plotExprHeatmap(sce_tcells, features = "type", by = "cluster_id", 
                k = "meta12", k_pal = colors_tcells, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()

