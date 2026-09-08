############################## README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: To normalize and scale data, compute variable features, 
# regress out e.g. nCount, 
# compute cell cycle scores,
# and run PCA.

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(readxl)

set.seed(100)

source("workflow/functions.R")

############################## Load tester variables
# rds_input <- "results/08_merge_samples/megaminus/merged_megaminus_DRscrub.rds"
# variables2regress <- c("nCount_RNA", "S.Score", "G2M.Score")
# ndims2check <- 50
# reduction <- "pca"

############################## Load Snakemake variables
rds_input <- snakemake@input[["rds"]]
variables2regress <- snakemake@params$variables2regress
ndims2check <- 50
reduction <- "pca"

plot7 <- snakemake@output[["plot7"]]
plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]
rds_output <- snakemake@output[["rds"]]

############################## Load data
sf <- readRDS(rds_input)

############################## Clean up object
sf <- JoinLayers(sf)
sf
sf <- DietSeurat(sf, dimreducs = NULL, layers = "counts")
sf
sf[["RNA"]] <- split(sf[["RNA"]], f = sf$orig.ident)

############################## Normalize data
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

############################## Calculate Cell Cycle Score
# Needs to be done before data is scaled, because the cell cycle scores are
# to be regressed 

# Read in compiled list of cell cycle genes
cellcycle_genes <- read_xlsx("miscellaneous/cellcycle_genes.xlsx")

# Isolate S-phase and G2/M phase genes
s.genes <- cellcycle_genes$genename_Hs_PwNameFormat[cellcycle_genes$phase == "S"]
g2m.genes <- cellcycle_genes$genename_Hs_PwNameFormat[cellcycle_genes$phase == "G2M"]

# Run Cell Cycle scoring. Runs on each layer separately
sflognorm <- CellCycleScoring(sflognorm, s.features = s.genes, g2m.features = g2m.genes)

# Visualize the Cell Cycle Score across samples
VlnPlot(sflognorm, c("S.Score", "G2M.Score"), group.by = "orig.ident")
ggsave(filename = plot7, width = 10)

############################## Check if variables2regress are available
print(variables2regress)

# If the variable is NULL, then set the value to NULL
if ("NULL" %in% variables2regress) {
  variables2regress <- NULL
} else {
  # If the variable is not NULL, check that it's present in the object
  if (!all(variables2regress %in% colnames(sflognorm@meta.data))) {
    stop(paste("The variable(s) to regress is not present in the object:",
               paste0(variables2regress[!variables2regress %in% colnames(sflognorm@meta.data)], collapse = ", ")))
  }
}

############################## Scale data

all.genes <- rownames(sflognorm)
sflognorm <- ScaleData(sflognorm, vars.to.regress = variables2regress)

############################## Dimensionality reduction
sflognorm <- RunPCA(sflognorm,
                    features = VariableFeatures(object = sflognorm),
                    seed.use = 100)

############################## Plotting
# PCA plot
DimPlot(sflognorm, reduction = "pca", group.by = "orig.ident", shuffle = TRUE)
ggsave(filename = plot1)

## nCount & nFeature overlaid on PCA plot
FeaturePlot(sflognorm,
       features = c("nCount_RNA", "nFeature_RNA"),
       reduction = "pca")
ggsave(filename = plot2, width = 10)

# Save cumulative PCs contribution
plot_elbow(sflognorm)
ggsave(filename = plot3, width = 15)

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

############################## Save the Log normalized data
saveRDS(sflognorm, rds_output)

############################## Close log file
sink(type = "message")
sink()
close(logfile)