# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aims:

# 1. For a given input dimensionality reduction (typically PCA, integrated.cca, etc.),
# get the UMAP dimensionality reduction for plotting purposes

# 2. Visualize known markers genes of interest & canonical markers
# on the UMAPs

# 3. Compute and plot clusters for a range of resolutions,
# and evaluate overlap with known markers, as well as cluster stability.

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(clustree)
library(gridExtra)
library(grid)
library(pals)
library(patchwork)

source("workflow/functions.R")

set.seed(100)
############################## Load tester variables
# ndims <- 25
# RunUMAP_input <- "integrated.cca"
# RunUMAP_output <- "umap.cca"
# cluster_prefix <- "cca"

############################## Load Snakemake variables
input_rds <- snakemake@input[["input_rds"]]
ndims <- as.numeric(snakemake@params$ndims) # User-defined ndims
RunUMAP_input <- snakemake@params$RunUMAP_input # User-defined dimensionality reduction
# Can be "pca", "integrated.cca", "integrated.rpca", "integrated.harmony", etc.
RunUMAP_output <- snakemake@params$RunUMAP_output
cluster_prefix <- snakemake@params$cluster_prefix

knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

txt <- snakemake@output[["txt"]]
plot0 <- snakemake@output[["plot0"]]
plot1_1 <- snakemake@output[["plot1_1"]]
plot1_2 <- snakemake@output[["plot1_2"]]
plot2_1 <- snakemake@output[["plot2_1"]]
plot2_2 <- snakemake@output[["plot2_2"]]
plot2_3 <- snakemake@output[["plot2_3"]]

plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
output_table <- snakemake@output[["output_table"]]
output_rds <- snakemake@output[["output_rds"]]

############################### Load data
sfi <- readRDS(input_rds)

# For clean feature plots
sfi <- JoinLayers(sfi)

############################### Check snakemake variables
RunUMAP_input %in% names(sfi@reductions)

if (!RunUMAP_input %in% names(sfi@reductions)) {
  stop(paste0(
    RunUMAP_input,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(sfi@reductions), collapse = ", ")
  ))
}

############################### Ndims and reductions
# Record the ndims and reduction used in an accompanying text file
writeLines(c(paste0("ndims: ", ndims),
             paste0("reduction input to RunUMAP:", RunUMAP_input),
             paste0("reduction output from RunUMAP:", RunUMAP_output),
             paste0("clusterings at different resolutions have the prefixes:",
                    cluster_prefix)),
           con = txt)

############################### UMAP dimensionality reduction

sfi <- RunUMAP(sfi,
               dims = 1:ndims,
               reduction = RunUMAP_input,
               reduction.name = RunUMAP_output)

############################### UMAP dimensionality reduction
# Split by sample
png(plot2_1)
DimPlot(sfi, reduction = RunUMAP_output)
dev.off()

# Split by batch
png(plot2_2)
DimPlot(sfi,
        reduction = RunUMAP_output,
        group.by = "batch",
        cols = c("purple", "darkgreen"))
dev.off()

# Split by condition
png(plot2_3)
DimPlot(sfi, reduction = RunUMAP_output, group.by = "condition")
dev.off()

################################ Basic feature plots
# Specify genes of interest (known marker genes)
goi <- read.csv(knownmarker_file, sep = ";")

goi <- goi[!is.na(goi[[knownmarker_id_col]]) & goi[[knownmarker_id_col]] != "", ]
goi <- goi[!duplicated(goi[[knownmarker_id_col]]), ]

names_goi <- goi[[knownmarker_id_col]]
goi2plot <- names_goi[names_goi %in% rownames(sfi)]

if (length(goi2plot) == 0) {
  warning("No marker genes from the marker file are present in the Seurat object. Please check your marker file and Seurat object rownames.")
  stop("No marker genes found in Seurat object. Exiting.")
}


# Split list of genes per page chunks
chunk_size <- 6
gene_chunks <- split(goi2plot, ceiling(seq_along(goi2plot) / chunk_size))

# Open PDF in landscape A4 size (11.7 x 8.3 inches)
pdf(plot0, width = 8.3, height = 11.7)

for (genes in gene_chunks) {
  p <- FeaturePlot(
    sfi,
    features = genes,
    reduction = RunUMAP_output,
    cols = c("grey", "darkred")
  )
  print(p)  # Each print = new page
}

dev.off()

############################### Test clustering resolutions
# Find neighbors using PCA
sfi <- FindNeighbors(sfi, reduction = RunUMAP_input, dims = 1:ndims)

# No. of clusters per resolution: calculate
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
num_clusters_per_resolution <- vector(length = length(resolutions))

# Desired order of images
new_order <- gsub("-", ".", levels(sfi$orig.ident))
sfi@images <- sfi@images[new_order]

# Algorithm to choose:
# 1 = original Louvain algorithm;
# 2 = Louvain algorithm with multilevel refinement;
# 3 = SLM algorithm;
# 4 = Leiden algorithm

umaps_clusters <- list()
spatial_clusters <- list()

for (i in seq_along(resolutions)) {

  # Find clusters at a given resolution
  sfi <- FindClusters(sfi, resolution = resolutions[i], algorithm = 4)

  # Record how many clusters came from this resolution
  num_clusters_per_resolution[i] <- length(unique(sfi$seurat_clusters))

  # adds 0 padding to make single-digit number like 9 into 09
  # for the cluster ids (per cell)
  ids <- as.character(sprintf("%02.0f", as.integer(sfi$seurat_clusters)))

  # turn the cluster ids into a factor, sorted in ascending order
  ids <- factor(ids, levels = sort(unique(ids)))

  # the vector of cluster ids should have the names of corresponding cell barcodes
  names(ids) <- colnames(sfi)

  # name of column to hold clustering results
  newcolname <- paste0(cluster_prefix, "_resolution_", as.character(resolutions[i]))

  # Add the padded cluster ids as a metadata column
  sfi[[newcolname]] <- ids

  # High-contrast color palette
  selected_cols <- pals::glasbey(30)
  names(selected_cols) <- Idents(sfi) %>% levels()

  # Make UMAP
  p <- DimPlot(sfi, label = TRUE, reduction = RunUMAP_output, cols = selected_cols) +
    ggtitle(newcolname) +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot all images with clustering info
  p_spatial <- plotall(sfi,
               pt.size = 3,
               selected_cols = selected_cols)

  # put all images on the same page
  combined_spatial <- wrap_plots(p_spatial, ncol = 3, guides = "collect")

  # Store the UMAP and Spatial plot
  umaps_clusters[[i]] <- p
  spatial_clusters[[i]] <- combined_spatial
}

# Save UMAPs of clusters
pdf(plot1_1)
umaps_clusters
dev.off()

# Save spatial plots of clusters
pdf(plot1_2, width = 8.3, height = 11.7)
spatial_clusters
dev.off()

############################### Save clustering evaluation plots
# No. of clusters per resolution: save table
resolutions_clusters <- data.frame(resolutions = resolutions,
                                   no.of.clusters = num_clusters_per_resolution)
write_xlsx(x = resolutions_clusters, path = output_table)

# No. of clusters per resolution: plot
p <- resolutions_clusters %>%
  ggplot(aes(x = resolutions, y = no.of.clusters)) +
  geom_smooth(se = FALSE, linetype = 2, colour = "grey20", linewidth = 0.4) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = round(resolutions, digits = 1)) +
  scale_y_continuous(breaks = unique(resolutions_clusters$no.of.clusters)) +
  labs(x = "Resolution", y = "Number of clusters")
p

ggsave(plot3, p, height = 4, width = 6)

## Cluster stability across resolutions
pdf(plot4, width = 14, height = 12)
ct <- clustree(sfi,
               prefix = paste0(cluster_prefix, "_resolution_"),
               layout = "sugiyama",
               exprs = "data") +
  guides(edge_colour = "none", edge_alpha = "none") +
  labs(y = NULL, x = NULL) +
  theme(legend.position = "right", legend.direction = "vertical")

# Create table grob from data frame for annotation
legend_tbl <- tableGrob(resolutions_clusters, rows = NULL)

# Arrange plot and table side by side
grid.arrange(ct, legend_tbl, ncol = 2, widths = c(4, 1))

ct_stability <- clustree(sfi,
                         prefix = paste0(cluster_prefix, "_resolution_"),
                         layout = "sugiyama",
                         exprs = "data",
                         node_colour = "sc3_stability") +
  guides(edge_colour = "none", edge_alpha = "none") +
  labs(y = NULL, x = NULL) +
  theme(legend.position = "right", legend.direction = "vertical")

ct_stability

dev.off()

################################ Save object
# Save the dimensionality reduced data
saveRDS(sfi, output_rds)

############################### Close log file
sink(type = "message")
sink()
close(logfile)