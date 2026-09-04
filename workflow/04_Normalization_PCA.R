############################## README
# Author: Inika Prasad
# Aim: To do normalization and dimensionality reduction

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

###############################  Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

set.seed(100)

source("workflow/functions.R")

############################### Tester variables
# ndims2check <- 30
# reduction <- "pca"

############################### Load Snakemake variables
input_rds <- snakemake@input[["input_rds"]]

ndims2check <- snakemake@params[["ndims2check"]]
reduction <- snakemake@params[["reduction"]]
variables2regress <- snakemake@params$variables2regress

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]
output_rds <- snakemake@output[["output_rds"]]

############################### Load data
sf <- readRDS(input_rds)

# If the variable is NULL, then set the value to NULL
if (variables2regress == "NULL") {
  variables2regress <- NULL
} else {
  # If the variable is not NULL, check that it's present in the object
  if (!(variables2regress %in% colnames(sf@meta.data))) {
    stop("The variable to regress is not present in the object")
  }
}

################################# Normalize data
sflognorm <- sf %>%
  NormalizeData() %>%
  FindVariableFeatures() %>% # ELEPHANT
  ScaleData()

# Identify the 10 most highly variable genes (probably a consensus set)
top10 <- head(VariableFeatures(sflognorm), 10)
message("These are the top 10 most highly variable features")
top10

# Plot top Variable Features
png(plot0, width = 900)
plot_VF <- VariableFeaturePlot(sflognorm)
plot_VF_labelled <- LabelPoints(plot = plot_VF, points = top10, repel = TRUE)
plot_VF + plot_VF_labelled
dev.off()

## Dimensionality reduction
sflognorm <- RunPCA(sflognorm,
                    # features = VariableFeatures(object = sflognorm),
                    seed.use = 100)

## PCA plot
DimPlot(sflognorm, reduction = "pca", group.by = "orig.ident", shuffle = TRUE)
ggsave(filename = plot1)

################################# Make & sabe plots

# nCount & nFeature overlaid on PCA plot
FeaturePlot(sflognorm,
       features = c("nCount_Spatial", "nFeature_Spatial"),
       reduction = "pca")
ggsave(filename = plot2,
       width = 10)

# Save cumulative PCs contribution
plot_elbow(sflognorm)
ggsave(filename = plot3,
       plot = plot_elbow(sflognorm),
       width = 15)

# Save individual PC contribution
ElbowPlot(sflognorm)
ggsave(filename = plot4,
       plot = ElbowPlot(sflognorm, ndims = ndims2check),
       width = 15,
       bg = "white")

# Visualize cells in each PC
pdf(plot5, height = ndims2check * 1.6, width = 15)
DimHeatmap(sflognorm, dims = 1:ndims2check, cells = 500, balanced = TRUE) +
  theme(plot.margin = margin(5.5, 20, 5.5, 10))
dev.off()

# Visualize genes & loadings in each PC
pdf(plot6, height = ndims2check * 1.6, width = 15)
VizDimLoadings(sflognorm, dims = 1:ndims2check, reduction = reduction) +
  theme(plot.margin = margin(5.5, 20, 5.5, 10))
dev.off()

################################# Save object
# Save the Log normalized data
saveRDS(sflognorm, output_rds)

############################### Close log file
sink(type = "message")
sink()
close(logfile)