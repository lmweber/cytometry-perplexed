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


# -------------------
# Dimension reduction
# -------------------

n_dims <- 2


# run Rtsne with different seeds

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


# run UMAP with different seeds

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

colors_clus <- c(rainbow(12))

# function for tSNE plots
plot_tSNE_clus <- function(dims_tSNE, labels, colors) {
  
  d_plot <- cbind(as.data.frame(dims_tSNE), cluster = as.factor(labels))
  
  ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

# seed 1
plot_tSNE_clus(dims_Rtsne_1, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed1.png", width = 4.5, height = 3.5)

# seed 2
plot_tSNE_clus(dims_Rtsne_2, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed2.png", width = 4.5, height = 3.5)

# seed 3
plot_tSNE_clus(dims_Rtsne_3, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed3.png", width = 4.5, height = 3.5)


# function for UMAP plots
plot_UMAP_clus <- function(dims_umap, labels, colors) {
  
  d_plot <- cbind(as.data.frame(dims_umap), cluster = as.factor(labels))
  
  ggplot(d_plot, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

# seed 1
plot_UMAP_clus(dims_umap_1, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed1.png", width = 4.5, height = 3.5)

# seed 2
plot_UMAP_clus(dims_umap_2, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed2.png", width = 4.5, height = 3.5)

# seed 3
plot_UMAP_clus(dims_umap_3, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed3.png", width = 4.5, height = 3.5)


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
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

# seed 1
plot_tSNE_truth(dims_Rtsne_1, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed1.png", width = 6, height = 3.5)

# seed 2
plot_tSNE_truth(dims_Rtsne_2, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed2.png", width = 6, height = 3.5)

# seed 3
plot_tSNE_truth(dims_Rtsne_3, pop_ids_sub, colors_truth)
ggsave("../plots/truth_tSNE_seed3.png", width = 6, height = 3.5)


# function for UMAP plots
plot_UMAP_truth <- function(dims_umap, pop_ids, colors) {
  
  d_plot <- cbind(as.data.frame(dims_umap), population = as.factor(pop_ids))
  
  ggplot(d_plot, aes(x = UMAP_1, y = UMAP_2, color = population)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = colors) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
}

# seed 1
plot_UMAP_truth(dims_umap_1, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed1.png", width = 6, height = 3.5)

# seed 2
plot_UMAP_truth(dims_umap_2, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed2.png", width = 6, height = 3.5)

# seed 3
plot_UMAP_truth(dims_umap_3, pop_ids_sub, colors_truth)
ggsave("../plots/truth_UMAP_seed3.png", width = 6, height = 3.5)


# -------------------------
# Clustering within T-cells
# -------------------------

# identify T-cells from FlowSOM clustering by inspecting UMAP and heatmap

# T cells are cluster 1
table(labels_fs)
# number and proportion of cells identified as T cells
table(labels_fs == 1)
mean(labels_fs == 1)

d_sub_tcells <- d_sub[labels_fs == 1, ]
dim(d_sub_tcells)

# perform second round of clustering on T cells only
d_flowsom_tcells <- flowFrame(d_sub_tcells)
set.seed(123)
out_fs_tcells <- ReadInput(d_flowsom_tcells, transform = FALSE, scale = FALSE)
out_fs_tcells <- BuildSOM(out_fs_tcells)
out_fs_tcells <- BuildMST(out_fs_tcells)
labels_pre_meta_tcells <- out_fs_tcells$map$mapping[, 1]
# number of meta-clusters
k <- 6
seed <- 123
out_fs_tcells <- metaClustering_consensus(out_fs_tcells$map$codes, k = k, seed = seed)
labels_fs_tcells <- out_fs_tcells[labels_pre_meta_tcells]
# summary of clusters
table(labels_fs_tcells)


# plot T cell subsets: tSNE

dims_Rtsne_1_tcells <- dims_Rtsne_1[labels_fs == 1, ]
stopifnot(nrow(dims_Rtsne_1_tcells) == length(labels_fs_tcells))

colors_clus <- rainbow(6)

plot_tSNE_clus(dims_Rtsne_1_tcells, labels_fs_tcells, colors_clus)
ggsave("../plots/FlowSOM_tSNE_seed1_Tcells.png", width = 4.5, height = 3.5)


# plot T cell subsets: UMAP

dims_umap_1_tcells <- dims_umap_1[labels_fs == 1, ]
stopifnot(nrow(dims_umap_1_tcells) == length(labels_fs_tcells))

colors_clus <- rainbow(6)

plot_UMAP_clus(dims_umap_1_tcells, labels_fs_tcells, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed1_Tcells.png", width = 4.5, height = 3.5)

