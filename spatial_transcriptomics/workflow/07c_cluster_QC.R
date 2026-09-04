#################################### README
# Author: Inika Prasad
# Created on: 2025-12-01
# Aim: EValuate the cluster QC metrics and across samples

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

################################ Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(writexl)
library(readxl)
library(patchwork)

library(ggpointdensity)
library(ggnewscale)

set.seed(100)

source("workflow/functions.R")

################################# Load tester variables
# Load data
input_rds <- "../results/06_dimred_clustering/allsections/allsections_dimred.rds"
myclusters <- "cca_resolution_1.2"
myreduction <- "umap.cca"

plot5 <- "condition_props_percluster.png"

################################# Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot1_1 <- snakemake@output[["plot1_1"]]
plot1_2 <- snakemake@output[["plot1_2"]]
plot4 <- snakemake@output[["plot4"]]
table1 <- snakemake@output[["table1"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]

myreduction <- snakemake@params$myreduction # UMAP to be used for plots
myclusters <- snakemake@params$myclusters # Clustering that is being QC'd
knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

################################# Load dataset
obj <- readRDS(input_rds)

################################# Check that inputs are valid
# Check that the chosen reduction and clustering exist,
# Abort with a message if they don't

if (!myreduction %in% names(obj@reductions)) {
  stop(paste0(
    myreduction,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(obj@reductions), collapse = ", ")
  ))
}

if (!myclusters %in% colnames(obj@meta.data)) {
  stop(paste0(
    myclusters,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  ))
}

if (!file.exists(knownmarker_file)) {
  stop(paste0(
    knownmarker_file, " does not exist. Please check path."
  ))
}

################################# Dotplot known marker genes
# Specify genes of interest (known marker genes)
goi <- read.csv(knownmarker_file, sep = ";")
print(goi)

# Clean up list
goi <- goi[goi[[knownmarker_id_col]] != "",]
goi <- goi[!duplicated(goi[[knownmarker_id_col]]), ]

# Create list of named lists
celltype_lists <- goi %>%
  group_by(Cell.Type) %>%
  summarise(gene = list(!!sym(knownmarker_id_col)), .groups = "drop") %>%
  { setNames(.$gene, .$Cell.Type) }

# Reorder list
categories_order <- unique(goi$Cell.Type)

celltype_lists <- setNames(
  lapply(categories_order, function(ct) goi[[knownmarker_id_col]][goi$Cell.Type == ct]),
  categories_order
)

pdf(plot0, width = 30, height = 8.3)

# Scaled plot
DotPlot(obj,
        group.by = myclusters,
        features = celltype_lists,
        cols = "RdBu",
        scale = TRUE) +
  RotatedAxis()

# Unscaled plot
DotPlot(obj,
        group.by = myclusters,
        features = celltype_lists,
        cols = "RdBu",
        scale = FALSE) +
  RotatedAxis()

dev.off()

################################## QC plots

# Plot sample of origin
p1 <- DimPlot(obj,
              reduction = myreduction,
              group.by = "orig.ident",
              label.size = 2) +
  ggtitle("Samples")
p1
ggsave(plot1)

# Plot clusters
p2 <- DimPlot(obj,
              reduction = myreduction,
              group.by = myclusters,
              label = TRUE) +
  ggtitle(myclusters)
p2
ggsave(plot2, width = 8)

# Split into conditions, monochrome
p1_1 <- DimPlot(obj,
                reduction = myreduction,
                split.by = "condition",
                group.by = "condition") +
  ggtitle("Samples") +
  theme(legend.position = "bottom")
p1_1
ggsave(plot1_1, width = 13)

# Split into samples, color by clusters
p1_2 <- DimPlot(obj,
                reduction = myreduction,
                split.by = "condition",
                group.by = myclusters,
                label = TRUE) +
  ggtitle(myreduction)
p1_2
ggsave(plot1_2, width = 18)

# Plot nCount and nFeature
p4 <- FeaturePlot(obj,
                  reduction = myreduction,
                  features = c("nCount_Spatial", "nFeature_Spatial"),
                  label.size = 2)
p4
ggsave(plot4, width = 16)

############################################ QC: no. of cells across clusters & samples
# No. of cells per sample per cluster
df <- table(obj@meta.data[["condition"]], obj@meta.data[[myclusters]])

# Save no. of cells per sample per cluster
write_xlsx(x = df %>% as.data.frame.matrix(), path = table1)

############################################ QC Plot: sample proportions across clusters
# all clusters in a sample add up to 100% (account for different no. of cells per sample)
propmat <- prop.table(df, margin = 1)

# proportions per cluster add up to 100% (account for cluster size)
propmat <- prop.table(propmat, margin = 2) * 100

df2plot <- propmat %>% as.data.frame()
colnames(df2plot) <- c("Sample", "Cluster", "Proportion")

# Get numbers and proportions of cells per cluster
print(df)
print(df2plot)

# hierachical clustering for heatmap
dist_cols <- dist(t(propmat)) # compute distance
hc_cols <- hclust(dist_cols) # cluster columns based on distance
category_order <- hc_cols$labels[hc_cols$order] #reorder labels based on clustering
df2plot$Cluster <- factor(df2plot$Cluster, levels = category_order) #reorder Cluster levels

# Heatmap of proportions
rd_bu_colors <- rev(brewer.pal(5, "RdBu")) # choosing 3 makes the heatmap paler
ggplot(df2plot, aes(x = Cluster, y = Sample, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = rd_bu_colors,
    values = scales::rescale(c(0, 33, 100)),  # midpoint at 0.25
    limits = c(0, 100),
    name = "Proportion (%)"
  ) +
  theme_minimal() +
  xlab(myclusters) +
  ggtitle("Sample proportions per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

nclust <- df2plot$Cluster %>% unique() %>% length()
ggsave(plot5, width = 6 + (nclust * 0.14), height = 2, bg = "white")

############################################ QC Plot: doublets per cluster
# proportion predicted doublets per cluster
df <- table(obj@meta.data[["Doublet_Confidence"]], obj@meta.data[[myclusters]])
propmat <- prop.table(df, margin = 2) * 100
df2plot <- propmat %>% as.data.frame()
colnames(df2plot) <- c("Doublet_Confidence", "Cluster", "Proportion")

# hierachical clustering for heatmap
dist_cols <- dist(t(propmat)) # compute distance
hc_cols <- hclust(dist_cols) # cluster columns based on distance
category_order <- hc_cols$labels[hc_cols$order] #reorder labels based on clustering
df2plot$Cluster <- factor(df2plot$Cluster, levels = category_order) #reorder Cluster levels

# Heatmap of proportions
rd_bu_colors <- rev(brewer.pal(5, "RdBu"))
ggplot(df2plot, aes(x = Cluster, y = Doublet_Confidence, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = rd_bu_colors,
    values = scales::rescale(c(0, 50, 100)),
    limits = c(0, 100),
    name = "Proportion (%)"
  ) +
  theme_minimal() +
  xlab(myclusters) +
  ggtitle("Doublet proportions per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

nclust <- df2plot$Cluster %>% unique() %>% length()
ggsave(plot6, width = 6 + (nclust * 0.14), height = 2, bg = "white")

#################################### Close log file
sink(type = "message")
sink()
close(logfile)