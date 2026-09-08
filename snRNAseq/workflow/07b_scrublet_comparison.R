############################# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Created on: 29 September, 2025
# Aim: To evaluate Scrublet results (doublet score and classification),
# Side by side with DoubletFinder scores

# If one were using the Scrublet from scratch, one might
# remove all the cells designated as doublets before running any filtering,
# variable gene selection, or even PCA & UMAP...
# But we will do that in the merged datasets anyways

############################# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################# Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(gridExtra)
library(grid)

source("workflow/functions.R")

############################ Set seed for reproducibility
set.seed(100)

############################ Load snakemake variables
rds_input <- snakemake@input[["rds_input"]]
csv_input <- snakemake@input[["csv_input"]]

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
plot5 <- snakemake@output[["plot5"]]
txt1 <- snakemake@output[["txt1"]]
rds_output <- snakemake@output[["rds_output"]]

chosen_resolution <- snakemake@params[["chosen_resolution"]]

############################# Load data
sfi <- readRDS(rds_input)
scrub_info <- read.csv(csv_input, sep = ";")

sample_name <- sfi$orig.ident %>% unique()

############################# Add Scrublet info to Seurat object
# Preview Scrublet info
print("Preview Scrublet info. Snippet of scrublet info to incorporate into the Seurat object:")
head(scrub_info)

# Check dimensions of the Scrublet info
print("No. of barcodes for which we have Scrublet info:")
nrow(scrub_info)

print("No. of barcodes in the Seurat object after filtering:")
length(colnames(sfi))

# Check that all the Seurat object barcodes are in the metadata
if (all(colnames(sfi) %in% scrub_info$barcode)) {
    print("Every barcode has a Scrublet prediction: Yes")
} else {
    print("Some barcodes in the Seurat object do not have a Scrublet prediction. 
           Ensure that the csv_input
           file with Scrublet prediction matches the Seurat object")
}

# Check if Scrublet predictions were made for more cells than we kept in the Seurat object
# Since I tried to make the Scrublet run on a count matrix with the same filtering
if (nrow(scrub_info) == length(colnames(sfi))){
    print("Scrublet predictions are made for the same number of cells. Excellent.")
} else {
    print("Scrublet made predictions for a different number of cells. 
           Check your Scrublet run to ensure you cleaned the count matrix in the same way")
}

# Add a barcode column to the Seurat metadata
sfi@meta.data$barcode <- colnames(sfi)

# Combine the existing metadata and the srublet info 
new_metadata <- left_join(sfi@meta.data, scrub_info, by = "barcode")

# # Inspect the new metadata
# head(new_metadata)

# Add the new metadata to the Seurat object
sfi <- AddMetaData(object = sfi, metadata = new_metadata)

############################# No. of doublets predicted

# How many doublets predicted in each group?
print("No. of doublets predicted by DoubletFinder")
tbl_doubletfinder <- table(sfi$scrublet_DropletType) %>% as.data.frame()
tbl_doubletfinder

print("No. of doublets predicted by Scrublet")
tbl_scrublet <- table(sfi$Doublet_Confidence) %>% as.data.frame()
tbl_scrublet 

# Create table grob from data frame for annotation
tbl_doubletfinder <- tableGrob(tbl_doubletfinder, rows = NULL)
tbl_scrublet <- tableGrob(tbl_scrublet, rows = NULL)

png(plot0)
# Arrange plot and table side by side
grid.arrange(tbl_doubletfinder, tbl_scrublet, ncol = 2)
dev.off()

############################# Scrublet & DoubletFinder: UMAP
DimPlot(sfi, group.by = c("Doublet_Confidence", "scrublet_DropletType"))
ggsave(filename = plot1, bg = "white", width = 20)

############################# Scrublet & DF: nCount & nFeature

# Visualize nCount & nFeature per for Scrublet groups
nCount_vs_nFeature <- FeatureScatter(sfi,
                                     feature1 = "nCount_RNA",
                                     feature2 = "nFeature_RNA",
                                     group.by = c("Doublet_Confidence", "scrublet_DropletType"),
                                     shuffle = TRUE)

ggsave(filename = plot3, bg = "white", width = 16)

############################# ScrubletScores: UMAP & Violin

# Visualize doublet scores from Scrublet
FeaturePlot(sfi, "scrublet_Scores", cols = c("white", "darkred"))
ggsave(filename = plot2, bg = "white")

# # Visualize Scrublet scores per classification
png(filename = plot4)
VlnPlot(sfi, features = "scrublet_Scores", group.by = chosen_resolution, alpha = 0.2) +
    ggtitle(paste0(sample_name, " Scrublet Scores"),
            subtitle = chosen_resolution)
dev.off()

############################# Doublet proportions per cluster
pdf(plot5)
plot_cluster_proportions(seurat_obj = sfi,
                         group_var = "scrublet_DropletType",
                         cluster_var = chosen_resolution,
                         sample_name = NULL)

plot_cluster_proportions(seurat_obj = sfi,
                         group_var = "Doublet_Confidence",
                         cluster_var = chosen_resolution,
                         sample_name = NULL)
dev.off()

############################# Save text summary:

# Information about what we saved
writeLines(c(paste0("The selected resolution:: ", chosen_resolution), 
             paste0("Clusters in the selected resolution:", paste(unique(sfi[[chosen_resolution]]), collapse = ", ")),
             paste0("No. of barcodes for which we have Scrublet info:", nrow(scrub_info)),
             paste0("No. of barcodes in the Seurat object after filtering:", length(colnames(sfi)))
             ), 
        con = txt1)

############################# Save object
# # Save object with Scrublet metadata included
saveRDS(sfi, rds_output)

#################################### Close log file
sink(type = "message")
sink()
close(logfile)