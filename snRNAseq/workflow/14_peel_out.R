# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Date created: 2025-09-19
# Aim: To take as input a Seurat object and a metadata field
# And to split the Seurat based on groupings of the metadata
# E.g. metadata cca_resolution_1.4 metadata column
# should take clusters 1,3,5,7 and put them in a new object
# which is saved as cardiomyocytes.rds

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")

################################ Load Libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

################################# Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

plot1 <- snakemake@output[["plot1"]]
txt1 <- snakemake@output[["txt1"]]
filename <- snakemake@output[["rds"]]

myclusters <- snakemake@params$myclusters # metadata column
clusters2group <- snakemake@params$clusters2group
cluster_groupname <- snakemake@wildcards[["cluster_name"]]
################################# Load dataset

obj <- readRDS(input_rds)

################################# Check if snakemake variables make sense
print(clusters2group)
if (!myclusters %in% colnames(obj@meta.data)) {
  stop(paste0(
    myclusters,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  ))
}

available_clusters <- obj@meta.data[[myclusters]] %>% unique()

if (!all(clusters2group %in% available_clusters)) {
  missing_clusters <- clusters2group[!clusters2group %in% available_clusters]
  stop(paste0(
        "The following clusters are not present in the metadata column '", myclusters, "': ",
        paste(missing_clusters, collapse = ", "), "\n",
        "Available metadata/clusterings are: ",
        paste(available_clusters, collapse = ", ")
  ))
}

################################# Subset dataset accordingly
# Assign the right clustering
Idents(obj) <- myclusters

# Subset out the cells of interest
obj_subset <- subset(x = obj, idents = clusters2group)
Idents(obj_subset) <- myclusters

# Ensure the resulting object has the exact clusters you require
# Visualize the cluster after subsetting
png(plot1)
DimPlot(obj_subset, label = TRUE)
dev.off()

# What clusters remain in the subsetted object?
clustersinsubset <- obj_subset@meta.data[[myclusters]] %>%
  unique() %>%
  as.character()

# Information about what we saved
writeLines(c(paste0("clustering used: ", myclusters),
             paste0("name of the group: ", cluster_groupname),
             paste0("clusters intended to be grouped: ", paste(clusters2group, collapse = ", ")),
             paste0("clusters actually grouped: ", paste(clustersinsubset, collapse = ", ")),
             paste0("name of resilting file: ", filename)),
           con = txt1)

# Save the object
saveRDS(obj_subset, filename)

#################################### Close log file
sink(type = "message")
sink()
close(logfile)