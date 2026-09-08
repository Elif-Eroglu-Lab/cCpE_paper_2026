# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(clustree)
library(gridExtra)
library(grid)

set.seed(100)

# Load snakemake variables
sample <- snakemake@params[["sample"]]
rds <- snakemake@input[["rds"]]
ndims <- as.numeric(snakemake@params$ndims) # User-defined ndims
reduction <- snakemake@params$reduction # User-defined dimensionality reduction

knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

txt1 <- snakemake@output[["txt"]]
plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
output_table <- snakemake@output[["table"]]
output_rds <- snakemake@output[["rds"]]

## Load data
sfi <- readRDS(rds)

# Record the ndims and reduction used in an accompanying text file
writeLines(c(paste0("ndims: ", ndims), paste0("reduction used:", reduction)), con = txt1)

# UMAP dimensionality reduction
sfi <- RunUMAP(sfi, dims = 1:ndims, reduction.name = "umap_chosen")
png(plot2)
DimPlot(sfi, reduction = "umap_chosen")
dev.off()

# Plot basic feature plots
# Specify genes of interest (known marker genes)
goi <- read.csv(knownmarker_file, sep = ";")
print(goi)
V2_goi <- goi[[knownmarker_id_col]]
print(V2_goi)
V2_goi <- V2_goi[V2_goi != ""] %>% unique()
print(V2_goi)

pdf(plot0, width = 50, height = 3.2 * length(V2_goi))
p <- FeaturePlot(sfi,
                  features = V2_goi,
                  reduction = "umap_chosen",
                  raster = TRUE)
print(p)
dev.off()

# Find neighbors using PCA
sfi <- FindNeighbors(sfi, reduction = reduction, dims = 1:ndims)

# No. of clusters per resolution: calculate
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
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

  # Add the padded cluster ids as a metadata column
  sfi[[paste0("resolution_", as.character(resolutions[i]))]] <- ids

  title4plot <- paste0("resolution:", as.character(resolutions[i]))
  p <- DimPlot(sfi, label = TRUE, reduction = "umap_chosen") + 
  ggtitle(title4plot) +
  theme(plot.title = element_text(hjust = 0.5))

  print(p)
}

dev.off()

# No. of clusters per resolution: save table
resolutions_clusters <- data.frame(resolutions = resolutions, no.of.clusters = num_clusters_per_resolution)
write_xlsx(x = resolutions_clusters, path = output_table)

# No. of clusters per resolution: plot
p <- resolutions_clusters %>%
  ggplot(aes(x = resolutions, y = no.of.clusters)) +
  geom_smooth(se = FALSE, linetype = 2, colour = "grey20", linewidth = 0.4) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = round(resolutions, digits = 1)) +
  scale_y_continuous(breaks = unique(resolutions_clusters$no.of.clusters)) +
  labs(x = "Resolution", y = "Number of clusters", title = sample)
p

ggsave(plot3, p, height = 4, width = 6)

## Cluster stability across resolutions
pdf(plot4, width = 14, height = 12)
ct <- clustree(sfi, prefix = "resolution_", layout = "sugiyama", exprs = "data") +
  guides(edge_colour = "none", edge_alpha = "none") +
  labs(y = NULL, x = NULL) +
  theme(legend.position = "right", legend.direction = "vertical")

# Create table grob from data frame for annotation
legend_tbl <- tableGrob(resolutions_clusters, rows = NULL)

# Arrange plot and table side by side
grid.arrange(ct, legend_tbl, ncol = 2, widths = c(4, 1))

ct_stability <- clustree(sfi, prefix = "resolution_", layout = "sugiyama", exprs = "data", node_colour = "sc3_stability") +
  guides(edge_colour = "none", edge_alpha = "none") +
  labs(y = NULL, x = NULL) +
  theme(legend.position = "right", legend.direction = "vertical")

ct_stability

dev.off()

## Save the dimensionality reduced data
saveRDS(sfi, output_rds)

# Close log file
sink(type = "message")
sink()
close(logfile)