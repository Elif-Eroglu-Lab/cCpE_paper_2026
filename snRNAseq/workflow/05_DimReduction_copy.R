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

set.seed(100)

############################## Load Snakemake variables
# sample <- snakemake@params[["sample"]]
rds <- snakemake@input[["rds"]]
ndims <- as.numeric(snakemake@params$ndims) # User-defined ndims
RunUMAP_input <- snakemake@params$RunUMAP_input # User-defined dimensionality reduction
# Can be "pca", "integrated.cca", "integrated.rpca", "integrated.harmony", etc.
RunUMAP_output <- snakemake@params$RunUMAP_output
cluster_prefix <- snakemake@params$cluster_prefix

knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

txt1 <- snakemake@output[["txt"]]
plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot2_1 <- snakemake@output[["plot2_1"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
output_table <- snakemake@output[["table"]]
output_rds <- snakemake@output[["rds"]]

############################## Load data
sfi <- readRDS(rds)

############################## Place samples in an intuitive order
sfi$orig.ident <- factor(sfi$orig.ident,
                         levels = c("Early_active",
                                    "PM_naive",
                                    "14dpsham",
                                    "10dpci_veh",
                                    "14dpci",
                                    "10dpci_cCpE"))

Idents(sfi) <- sfi$orig.ident

############################## Check variables2regress
# Check if the RunUMAP input exists in the object
RunUMAP_input %in% names(sfi@reductions)

if (!RunUMAP_input %in% names(sfi@reductions)) {
  stop(paste0(
    RunUMAP_input,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(sfi@reductions), collapse = ", ")
  ))
}


# Record the ndims and reduction used in an accompanying text file
writeLines(c(paste0("ndims: ", ndims),
             paste0("reduction input to RunUMAP:", RunUMAP_input),
             paste0("reduction output from RunUMAP:", RunUMAP_output),
             paste0("clusterings at different resolutions have the prefixes:",
                    cluster_prefix)),
           con = txt1)

############################## UMAP dimensionality reduction
sfi <- RunUMAP(sfi, dims = 1:ndims, reduction = RunUMAP_input, reduction.name = RunUMAP_output)
png(plot2)
DimPlot(sfi, reduction = RunUMAP_output)
dev.off()

# Split into samples, monochrome
p2_1 <- DimPlot(sfi,
                reduction = RunUMAP_output,
                split.by = "orig.ident",
                group.by = "orig.ident") +
  ggtitle("Samples") +
  theme(legend.position = "bottom")
p2_1
ggsave(plot2_1, width = length(unique(sfi$orig.ident)) * 5)

############################## Plot basic feature plots
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
    ncol = 2,
    features = genes,
    reduction = RunUMAP_output,
    raster = TRUE
  )
  print(p)  # Each print = new page
}

dev.off()

############################### Try different clustering resolutions
# Find neighbors using PCA
sfi <- FindNeighbors(sfi, reduction = RunUMAP_input, dims = 1:ndims)

# No. of clusters per resolution: calculate
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
num_clusters_per_resolution <- vector(length = length(resolutions))

# Algorithm to choose:
# 1 = original Louvain algorithm;
# 2 = Louvain algorithm with multilevel refinement;
# 3 = SLM algorithm;
# 4 = Leiden algorithm

pdf(plot1)

for (i in 1:length(resolutions)) {

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

  p <- DimPlot(sfi, label = TRUE, reduction = RunUMAP_output) +
    ggtitle(newcolname) +
    theme(plot.title = element_text(hjust = 0.5))

  print(p)
}

dev.off()

############################## Summary of clusters
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

############################## Save object
# Save data
saveRDS(sfi, output_rds)

############################## Close log file
sink(type = "message")
sink()
close(logfile)