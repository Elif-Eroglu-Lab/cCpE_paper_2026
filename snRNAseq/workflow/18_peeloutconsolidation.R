############################## README

# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: To concolidate the clustering labels from the peelouts 
# into one gigantic UMAP
# Date created: 2025-12-05

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

##################################### Load libraries
library(dplyr)
library(Seurat)
library(readxl)

##################################### Define the broad clustering dataset
fulldatafile <- "results/12_integrated_DimReduction/10dpci_ctrl/DRscrub/10dpci_ctrl_DRscrub_dimred.rds"
fulldatares <- "harmony_resolution_0.6"

consolidated_rds <- "results/testerdump/10dpci_ctrl_DRscrub_consolidated.rds"
#################################### Load the main dataset
fulldata <- readRDS(fulldatafile)

# Create a new cluster identity to hold main clusters & peelout clusters
fulldata$new_ident <- as.character(fulldata@meta.data[[fulldatares]])
fulldata$new_ident %>% unique()

# Visualize the clusters in single-nuclei
DimPlot(fulldata, group.by = "new_ident", label = TRUE)

###################################### Define the peelout datasets

peelout_datasets <- c("results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred.rds")

peelout_resolutions <- c("harmony_peeled_resolution_1")

peelout_names <- c("Epi")

peelout_reductions <- c("umap_peeled.harmony")

# MAKE SURE THE CLUSTERING MATCHES WHAT YOU HAVE IN THIS FILE!! MAYBE NAME THE SHEET BY RESOLUTION
peelout_namefiles <- c("results/16_peeled_integrated_DimRed/10dpci_ctrl/Epi_custom_clusternames.xlsx")


####################################### Consolidate the peelouts
Idents(fulldata) <- fulldatares
fulldata$new_ident <- as.character(Idents(fulldata))

# Load the peelout datasets
for (i in seq_along(peelout_datasets)) {

  peelout_data <- readRDS(peelout_datasets[i])
  peelout_res <- peelout_resolutions[i]
  peelout_name <- peelout_names[i]
  peelout_reduction <- peelout_reductions[i]
  peelout_mapping <- read_xlsx(peelout_namefiles[i])

  # Rename peelout_clusters as needed
  mapping_vec <- setNames(peelout_mapping$clustername, peelout_mapping$cluster)
  peelout_data@meta.data[[peelout_res]] <- dplyr::recode(peelout_data@meta.data[[peelout_res]], !!!mapping_vec)

  # Plot peelout clusters
  print(DimPlot(peelout_data, group.by = peelout_res, reduction = peelout_reduction, label = TRUE, pt.size = 1.5))

  # Merge the peelout clustering with the main clustering
  fulldata@meta.data[Cells(peelout_data), "new_ident"] <- as.character(peelout_data@meta.data[[peelout_res]])

}

# Visualize the consolidated clusters
DimPlot(fulldata, group.by = "new_ident", label = TRUE, repel = TRUE)

####################################### Save the consolidated dataset
saveRDS(fulldata, consolidated_rds)
