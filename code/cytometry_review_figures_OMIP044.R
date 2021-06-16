#######################################################
# Script to generate figures for cytometry review paper
# Lukas Weber, June 2021
#######################################################

library(flowCore)
library(SummarizedExperiment)
library(FlowSOM)
library(Rtsne)
library(uwot)
library(ggplot2)


# ---------
# Load data
# ---------

d_fs <- read.FCS("../../data/OMIP-044 full stain CD45 live events.fcs", 
                 transformation = FALSE, truncate_max_range = FALSE)

# convert to SummarizedExperiment

dim(d_fs)

markers <- unname(parameters(d_fs)$desc)

col_data <- DataFrame(
  channel_name = colnames(d_fs), 
  marker_name = markers, 
  marker_class = c(rep("none", 6), rep("type", 14), "none", rep("type", 13), "none")
)

row_data <- DataFrame(sample = rep("sample01", nrow(d_fs)))

d_se <- SummarizedExperiment(
  rowData = row_data, 
  colData = col_data, 
  assays = list(exprs = exprs(d_fs))
)

col_names <- col_data$marker_name
col_names[is.na(col_names)] <- "NA"
colnames(d_se) <- col_names


# -------------
# Preprocessing
# -------------

# select cell type marker columns
d_clus <- assay(d_se)[, colData(d_se)$marker_class == "type"]
dim(d_clus)
d_clus[1:6, 1:6]

# transform data
cofactor <- 250
d_clus <- asinh(d_clus / cofactor)
summary(d_clus)

# subsample for figures
n <- 100000
set.seed(123)
ix <- sample(seq_len(nrow(d_clus)), n)
d_sub <- d_clus[ix, ]
dim(d_sub)

# remove any near-duplicate rows (required by Rtsne)
dups <- duplicated(d_sub)
table(dups)
d_sub <- d_sub[!dups, ]
dim(d_sub)

ix_keep <- ix[!dups]


# -------------------
# Dimension reduction
# -------------------

n_dims <- 2


# skip tSNE since too slow with large number of cells

# run Rtsne with different random seeds

# set.seed(101)
# out_Rtsne_1 <- Rtsne(as.matrix(d_sub), dims = n_dims)
# dims_Rtsne_1 <- out_Rtsne_1$Y
# colnames(dims_Rtsne_1) <- c("tSNE_1", "tSNE_2")
# head(dims_Rtsne_1)
# 
# set.seed(102)
# out_Rtsne_2 <- Rtsne(as.matrix(d_sub), dims = n_dims)
# dims_Rtsne_2 <- out_Rtsne_2$Y
# colnames(dims_Rtsne_2) <- c("tSNE_1", "tSNE_2")
# head(dims_Rtsne_2)
# 
# set.seed(103)
# out_Rtsne_3 <- Rtsne(as.matrix(d_sub), dims = n_dims)
# dims_Rtsne_3 <- out_Rtsne_3$Y
# colnames(dims_Rtsne_3) <- c("tSNE_1", "tSNE_2")
# head(dims_Rtsne_3)


# run UMAP with different random seeds

set.seed(111)
dims_umap_1 <- umap(d_sub)
colnames(dims_umap_1) <- c("UMAP_1", "UMAP_2")
head(dims_umap_1)

set.seed(222)
dims_umap_2 <- umap(d_sub)
colnames(dims_umap_2) <- c("UMAP_1", "UMAP_2")
head(dims_umap_2)

set.seed(333)
dims_umap_3 <- umap(d_sub)
colnames(dims_umap_3) <- c("UMAP_1", "UMAP_2")
head(dims_umap_3)


# ----------
# Clustering
# ----------

# clustering using FlowSOM

d_flowsom <- flowFrame(d_sub)

set.seed(100)

out_fs <- ReadInput(d_flowsom, transform = FALSE, scale = FALSE)
out_fs <- BuildSOM(out_fs)
out_fs <- BuildMST(out_fs)

labels_pre_meta <- out_fs$map$mapping[, 1]

# number of meta-clusters
k <- 20

seed <- 123
out_fs <- metaClustering_consensus(out_fs$map$codes, k = k, seed = seed)

labels_fs <- out_fs[labels_pre_meta]

# summary of clusters
table(labels_fs)


# --------------------------
# Generate plots: clustering
# --------------------------

colors_clus <- unname(palette.colors(palette = "Polychrome 36"))

#library(RColorBrewer)
#colors_clus <- brewer.pal(12, "Paired")

#library(colorspace)
#set.seed(7)
#colors_clus <- sample(qualitative_hcl(12))
#colors_clus <- sample(qualitative_hcl(20))

#colors_clus <- rainbow(12)
#colors_clus <- brewer.pal(12, "Set3")

# # function for tSNE plots
# plot_tSNE_clus <- function(dims_tSNE, labels, colors, 
#                            size = 0.5, alpha = 0.1) {
#   
#   d_plot <- cbind(as.data.frame(dims_tSNE), cluster = as.factor(labels))
#   
#   ggplot(d_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
#     geom_point(size = size, alpha = alpha) + 
#     scale_color_manual(values = colors) + 
#     theme_bw() + 
#     theme(legend.key.size = unit(0.75, "lines"), 
#           panel.grid = element_blank(), 
#           axis.text = element_blank(), 
#           axis.ticks = element_blank()) + 
#     guides(color = guide_legend(override.aes = list(size = 1.5, alpha = 0.75)))
# }
# 
# # seed 1
# plot_tSNE_clus(dims_Rtsne_1, labels_fs, colors_clus)
# ggsave("../plots/FlowSOM_tSNE_seed1.png", width = 4.5, height = 3.75)
# ggsave("../plots/FlowSOM_tSNE_seed1.pdf", width = 4.5, height = 3.75)
# 
# # seed 2
# plot_tSNE_clus(dims_Rtsne_2, labels_fs, colors_clus)
# ggsave("../plots/FlowSOM_tSNE_seed2.png", width = 4.5, height = 3.75)
# ggsave("../plots/FlowSOM_tSNE_seed2.pdf", width = 4.5, height = 3.75)
# 
# # seed 3
# plot_tSNE_clus(dims_Rtsne_3, labels_fs, colors_clus)
# ggsave("../plots/FlowSOM_tSNE_seed3.png", width = 4.5, height = 3.75)
# ggsave("../plots/FlowSOM_tSNE_seed3.pdf", width = 4.5, height = 3.75)


# function for UMAP plots
plot_UMAP_clus <- function(dims_umap, labels, colors, 
                           size = 0.2, alpha = 0.1) {
  
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
ggsave("../plots/FlowSOM_UMAP_seed1.png", width = 4.5, height = 3.75)
ggsave("../plots/FlowSOM_UMAP_seed1.pdf", width = 4.5, height = 3.75)

# seed 2
plot_UMAP_clus(dims_umap_2, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed2.png", width = 4.5, height = 3.75)
ggsave("../plots/FlowSOM_UMAP_seed2.pdf", width = 4.5, height = 3.75)

# seed 3
plot_UMAP_clus(dims_umap_3, labels_fs, colors_clus)
ggsave("../plots/FlowSOM_UMAP_seed3.png", width = 4.5, height = 3.75)
ggsave("../plots/FlowSOM_UMAP_seed3.pdf", width = 4.5, height = 3.75)


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
  geom_point(size = 0.1, alpha = 0.1) + 
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
  geom_point(size = 0.1, alpha = 0.1) + 
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


# --------------------------
# Cluster expression heatmap
# --------------------------

# use heatmap plotting function from CATALYST

library(CATALYST)
library(scales)

# format into SingleCellExperiment for CATALYST
sce <- SingleCellExperiment(
  rowData = colData(d_se)[colData(d_se)$marker_class == "type", ], 
  colData = rowData(d_se)[ix_keep, , drop = FALSE], 
  assays = list(
    exprs = t(d_sub)
  )
)

colData(sce)$cluster_id <- labels_fs

# structure required by plotExprHeatmap
colData(sce)$sample_id <- "sample01"
metadata(sce)$cluster_codes <- data.frame(som100 = factor(1:k), meta = factor(1:k))

# order markers alphabetically
sce <- sce[order(rowData(sce)$marker_name), ]

pdf("../plots/heatmap_clusters_global.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce, features = "type", by = "cluster_id", scale = "first", 
                k = "meta", k_pal = colors_clus, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()

pdf("../plots/heatmap_clusters_channel.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce, features = "type", by = "cluster_id", scale = "last", 
                k = "meta", k_pal = colors_clus, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()


# -------------------------
# Clustering within T cells
# -------------------------

# identify T cells from FlowSOM clustering by inspecting UMAP and heatmap

# T cells are clusters 7, 9:11, 13:18
table(labels_fs)
# number and proportion of cells identified as T cells
table(labels_fs %in% c(7, 9:11, 13:18))
mean(labels_fs %in% c(7, 9:11, 13:18))

d_sub_tcells <- d_sub[labels_fs %in% c(7, 9:11, 13:18), ]
dim(d_sub_tcells)

# perform second round of clustering on T cells only
d_flowsom_tcells <- flowFrame(d_sub_tcells)
set.seed(111)
out_fs_tcells <- ReadInput(d_flowsom_tcells, transform = FALSE, scale = FALSE)
out_fs_tcells <- BuildSOM(out_fs_tcells)
out_fs_tcells <- BuildMST(out_fs_tcells)
labels_pre_meta_tcells <- out_fs_tcells$map$mapping[, 1]
# number of meta-clusters
#k <- 12
k <- 20
seed <- 123
out_fs_tcells <- metaClustering_consensus(out_fs_tcells$map$codes, k = k, seed = seed)
labels_fs_tcells <- out_fs_tcells[labels_pre_meta_tcells]
# summary of clusters
table(labels_fs_tcells)


# calculate UMAP dimension reduction: T cells

n_dims <- 2

set.seed(123)
dims_umap_tcells <- umap(d_sub_tcells)
colnames(dims_umap_tcells) <- c("UMAP_1", "UMAP_2")
stopifnot(nrow(dims_umap_tcells) == length(labels_fs_tcells))


# plot T cells subsets: UMAP

colors_tcells <- colors_clus

plot_UMAP_clus(dims_umap_tcells, labels_fs_tcells, colors_tcells)
ggsave("../plots/Tcells_FlowSOM_UMAP.png", width = 4.5, height = 3.75)
ggsave("../plots/Tcells_FlowSOM_UMAP.pdf", width = 4.5, height = 3.75)


# heatmap: T cells

ix_tcells <- labels_fs %in% c(7, 9:11, 13:18)
stopifnot(length(ix_tcells) == length(ix_keep))
stopifnot(nrow(d_sub_tcells) == nrow(rowData(d_se)[ix_keep, , drop = FALSE][ix_tcells, , drop = FALSE]))

# format into SingleCellExperiment for CATALYST
sce_tcells <- SingleCellExperiment(
  rowData = colData(d_se)[colData(d_se)$marker_class == "type", ], 
  colData = rowData(d_se)[ix_keep, , drop = FALSE][ix_tcells, , drop = FALSE], 
  assays = list(
    exprs = t(d_sub_tcells)
  )
)

colData(sce_tcells)$cluster_id <- labels_fs_tcells

# structure required by plotExprHeatmap
colData(sce_tcells)$sample_id <- "sample01"
metadata(sce_tcells)$cluster_codes <- data.frame(som100 = factor(1:k), meta = factor(1:k))

# order markers alphabetically
sce_tcells <- sce_tcells[order(rowData(sce_tcells)$marker_name), ]

pdf("../plots/heatmap_Tcells_global.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce_tcells, features = "type", by = "cluster_id", scale = "first", 
                k = "meta", k_pal = colors_tcells, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()

pdf("../plots/heatmap_Tcells_channel.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce_tcells, features = "type", by = "cluster_id", scale = "last", 
                k = "meta", k_pal = colors_tcells, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()


# -------------------------------
# Clustering within DCs + myeloid
# -------------------------------

# identify DCs + myeloid from FlowSOM clustering by inspecting UMAP and heatmap

# DCs + myeloid are clusters 1, 2, 3, 5, 6
table(labels_fs)
# number and proportion of cells identified as DCs + myeloid
table(labels_fs %in% c(1, 2, 3, 5, 6))
mean(labels_fs %in% c(1, 2, 3, 5, 6))

d_sub_dcs <- d_sub[labels_fs %in% c(1, 2, 3, 5, 6), ]
dim(d_sub_dcs)

# perform second round of clustering on DCs + myeloid only
d_flowsom_dcs <- flowFrame(d_sub_dcs)
set.seed(111)
out_fs_dcs <- ReadInput(d_flowsom_dcs, transform = FALSE, scale = FALSE)
out_fs_dcs <- BuildSOM(out_fs_dcs)
out_fs_dcs <- BuildMST(out_fs_dcs)
labels_pre_meta_dcs <- out_fs_dcs$map$mapping[, 1]
# number of meta-clusters
#k <- 12
k <- 20
seed <- 123
out_fs_dcs <- metaClustering_consensus(out_fs_dcs$map$codes, k = k, seed = seed)
labels_fs_dcs <- out_fs_dcs[labels_pre_meta_dcs]
# summary of clusters
table(labels_fs_dcs)


# calculate UMAP dimension reduction: DCs + myeloid

n_dims <- 2

set.seed(123)
dims_umap_dcs <- umap(d_sub_dcs)
colnames(dims_umap_dcs) <- c("UMAP_1", "UMAP_2")
stopifnot(nrow(dims_umap_dcs) == length(labels_fs_dcs))


# plot DCs + myeloid subsets: UMAP

colors_dcs <- colors_clus

plot_UMAP_clus(dims_umap_dcs, labels_fs_dcs, colors_dcs)
ggsave("../plots/DCs_FlowSOM_UMAP.png", width = 4.5, height = 3.75)
ggsave("../plots/DCs_FlowSOM_UMAP.pdf", width = 4.5, height = 3.75)


# heatmap: DCs + myeloid

ix_dcs <- labels_fs %in% c(1, 2, 3, 5, 6)
stopifnot(length(ix_dcs) == length(ix_keep))
stopifnot(nrow(d_sub_dcs) == nrow(rowData(d_se)[ix_keep, , drop = FALSE][ix_dcs, , drop = FALSE]))

# format into SingleCellExperiment for CATALYST
sce_dcs <- SingleCellExperiment(
  rowData = colData(d_se)[colData(d_se)$marker_class == "type", ], 
  colData = rowData(d_se)[ix_keep, , drop = FALSE][ix_dcs, , drop = FALSE], 
  assays = list(
    exprs = t(d_sub_dcs)
  )
)

colData(sce_dcs)$cluster_id <- labels_fs_dcs

# structure required by plotExprHeatmap
colData(sce_dcs)$sample_id <- "sample01"
metadata(sce_dcs)$cluster_codes <- data.frame(som100 = factor(1:k), meta = factor(1:k))

# order markers alphabetically
sce_dcs <- sce_dcs[order(rowData(sce_dcs)$marker_name), ]

pdf("../plots/heatmap_DCs_global.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce_dcs, features = "type", by = "cluster_id", scale = "first", 
                k = "meta", k_pal = colors_dcs, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()

pdf("../plots/heatmap_DCs_channel.pdf", width = 6.75, height = 4.5)
plotExprHeatmap(sce_dcs, features = "type", by = "cluster_id", scale = "last", 
                k = "meta", k_pal = colors_dcs, hm_pal = viridis_pal()(6), 
                row_clust = FALSE, col_clust = FALSE)
dev.off()

