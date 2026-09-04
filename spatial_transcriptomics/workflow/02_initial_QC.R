############################## README
# Author: Inika Prasad
# Aim: Compute baseline QC for samples

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(patchwork)

source("workflow/functions.R")

set.seed(100)

############################### Load Snakemake variables
input_rds <- snakemake@input[["input_rds"]]
sample <- snakemake@wildcards[["sample"]]

plot1 <- snakemake@output[["plot1"]]
plot1_1 <- snakemake@output[["plot1_1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
table <- snakemake@output[["table"]]
output_rds <- snakemake@output[["output_rds"]]

############################### Load data
obj <- readRDS(input_rds)
############################## QC plot & save
# Violin plots grouped by sample
VlnPlot(obj,
        features = c("nFeature_Spatial", "nCount_Spatial"),
        ncol = 2, alpha = 0.1,
        group.by = "sample")

ggsave(filename = plot1, bg = "transparent", width = 8)

# Violin plots, grouped by condition
VlnPlot(obj,
        features = c("nFeature_Spatial", "nCount_Spatial"),
        ncol = 2, alpha = 0.1,
        group.by = "condition")

ggsave(filename = plot1_1, bg = "transparent")

# Split by condition
## Plot nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(obj,
                                     feature1 = "nCount_Spatial",
                                     feature2 = "nFeature_Spatial",
                                     group.by = "sample",
                                     split.by = "condition",
                                     plot.cor = FALSE,
                                     shuffle = TRUE)
nCount_vs_nFeature
ggsave(filename = plot2, bg = "transparent", width = 7, units = "in")

# Visualize nCount & nFeature spatially
nFeatplots <- plotall(obj,
             features2plot = "nFeature_Spatial",
             crop = TRUE,
             pt.size= 3)

nCountplots <- plotall(obj,
             features2plot = "nCount_Spatial",
             crop = TRUE,
             pt.size = 3)

pdf(plot3, width = 8.3, height = 11.7)
wrap_plots(nFeatplots, ncol = 3, guides = "collect")
wrap_plots(nCountplots, ncol = 3, guides = "collect")
dev.off()

## How many cells express a given gene? Help decide min.cells threshold
su_countmat <- LayerData(obj, assay = "Spatial", layer = "counts")
gene_counts <- rowSums(su_countmat > 0)

png(plot4)
hist(gene_counts, breaks = 100, main = "Number of Spots Expressing Each Gene")
dev.off()

obj <- JoinLayers(obj)

## QC metrics per sample
su_summary <- obj[[]] %>%
  group_by(sample) %>%
  summarise(
    n_spots = n(),
    mean_nCount = mean(nCount_Spatial, na.rm = TRUE),
    mean_nFeatures = mean(nFeature_Spatial, na.rm = TRUE)
  )

## Save QC metrics
write_xlsx(x = su_summary, path = table)

############################### Close log file
sink(type = "message")
sink()
close(logfile)