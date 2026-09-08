# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

set.seed(100)

# Load snakemake variables
rds_input <- snakemake@input[["rds"]]
variables2regress <- snakemake@params$variables2regress
plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]
rds_output <- snakemake@output[["rds"]]
ndims2check <- 50
reduction = "pca"

## Define Elbow plot function
#' @title Plot elbow plot
#' @description Creates elbow plot for a dimensionality reduction
#' @param object The Seurat object for which you want the elbow plot
#' @param dimred The dimensionality reduction. E.g. "pca"
plot_elbow <- function(object, dimred = "pca") {
  pct <- object[[dimred]]@stdev / sum(object[[dimred]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90%
  # and % variation associated with the PC as less than 5
  # to identify when getting more PCs has diminishing returns.
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  # Minimum of the two calculation
  pcs <- min(co1, co2)
  message(paste0("elbow: ", pcs))

  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

  # Elbow plot to visualize
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_minimal()+
    theme(plot.background = element_rect(fill = "white", color = NA))

  return(p)
}

## Load data
sf <- readRDS(rds_input)

# If the variable is NULL, then set the value to NULL
if (variables2regress == "NULL") {
  variables2regress <- NULL
} else {
  # If the variable is not NULL, check that it's present in the object
  if (!(variables2regress %in% colnames(sf@meta.data))) {
    stop("The variable to regress is not present in the object")
  }
}

## Normalize data
sflognorm <- NormalizeData(sf)
sflognorm <- FindVariableFeatures(sflognorm,
                                  selection.method = "vst",
                                  nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sflognorm), 10)
message("These are the top 10 most highly variable features")
top10

# Plot top Variable Features
png(plot0, width = 900)
plot_VF <- VariableFeaturePlot(sflognorm)
plot_VF_labelled <- LabelPoints(plot = plot_VF, points = top10, repel = TRUE)
plot_VF + plot_VF_labelled
dev.off()

## Scale data
all.genes <- rownames(sflognorm)
sflognorm <- ScaleData(sflognorm, vars.to.regress = variables2regress)

## Dimensionality reduction
sflognorm <- RunPCA(sflognorm,
                    features = VariableFeatures(object = sflognorm),
                    seed.use = 100)

## PCA plot
ggsave(filename = plot1,
       plot = DimPlot(sflognorm, reduction = "pca", group.by = "orig.ident", shuffle = TRUE))

## nCount & nFeature overlaid on PCA plot
ggsave(filename = plot2,
       plot = FeaturePlot(sflognorm,
       features = c("nCount_RNA", "nFeature_RNA"),
       reduction = "pca"),
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

## Save the Log normalized data
saveRDS(sflognorm, rds_output)

# Close log file
sink(type = "message")
sink()
close(logfile)