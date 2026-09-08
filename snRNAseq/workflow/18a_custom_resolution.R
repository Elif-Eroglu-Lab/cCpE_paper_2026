############################ README ######################################
# This was run in conjunction with 12b_clustermarkers_4intdata.R
# To identify top markers of the epi peelout with a custom resolution
# i.e. after merging certain clusters.

########################### Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(tibble)
library(readxl)

set.seed(100)
source("workflow/functions.R")

############################ Load variables ###############################

# Run this instead of lines 25-35, then proceed with lines 35-84
{ # 10dpci_ctrl Epi Peelout
dataset <- "10dpci_ctrl_EpiPeelout"
rds_input <- "results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred.rds"
merged_object <- "results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"

table1 <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1_markers.xlsx" # Excel of all markers
plot0 <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1_umap.png" # UMAP showing clustering
plot1 <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1_marker_dotplot.pdf" # Dotplots

myreduction <- "umap_peeled.harmony" # UMAP to be used for plots
myclusters <- "harmony_peeled_resolution_1" # Clustering that is being QC'd
# grouping.var <- snakemake@params$grouping_var # metadata parameter for FindConservedMarkers
calc_conserved_markers <- FALSE
}


{ # Megaminus Epi Peelout
dataset <- "megaminus_EpiPeelout"
rds_input <- "results/16_peeled_integrated_DimRed/megaminus/megaminus_DRscrub_Epi_dimred.rds"
merged_object <- "results/16_peeled_integrated_DimRed/megaminus/megaminus_DRscrub_Epi_dimred_CUSTOMRES2.rds"

table1 <- "results/16b_clustermarkers_4peeledintdata/megaminus/Epi_custom2_res1.6_markers.xlsx" # Excel of all markers
plot0 <- "results/16b_clustermarkers_4peeledintdata/megaminus/Epi_custom2_res1.6_umap.png" # UMAP showing clustering
plot1 <- "results/16b_clustermarkers_4peeledintdata/megaminus/Epi_custom2_res1.6_marker_dotplot.pdf" # Dotplots

myreduction <- "umap_peeled.harmony" # UMAP to be used for plots
myclusters <- "harmony_peeled_resolution_1.6" # Clustering that is being QC'd
# grouping.var <- snakemake@params$grouping_var # metadata parameter for FindConservedMarkers
calc_conserved_markers <- FALSE
}

################################# Load dataset ################################
obj <- readRDS(rds_input)
obj <- JoinLayers(obj)

# Assign identities to cells based on the present cluster
Idents(obj) <- myclusters

# What levels does the clustering have?
clusters <- levels(obj)

# How many clusters do we have?
clusters %>% length()

DimPlot(obj, label = TRUE)

############################### Merge clusters ###############################
if (dataset == "10dpci_ctrl_EpiPeelout") {
  obj@meta.data[["merged_clusters"]] <- recode(as.character(obj@meta.data[[myclusters]]),
                                               `01` = "EPDC_A",
                                               `05` = "EPDC_A",
                                               `06` = "EPDC_A",
                                               `08` = "EPDC_A")
  print("Merged clusters in 10dpci_ctrl_EpiPeelout")
}

if (dataset == "megaminus_EpiPeelout") {
  obj@meta.data[["merged_clusters"]] <- recode(as.character(obj@meta.data[[myclusters]]),
                                               `01` = "EPDC_B",
                                               `02` = "EPDC_B",
                                               `10` = "EPDC_B",
                                               `04` = "SMC",
                                               `06` = "SMC",
                                               `13` = "SMC")
  print("Merged clusters in megaminus_EpiPeelout")
}

############################### QC the merging #############################

# Visualize the merge
DimPlot(obj, group.by = "merged_clusters", label = TRUE)

# Trying this out to see if it helps the NA-related bug in the DGE step
obj$orig.ident <- factor(obj$orig.ident, levels = c("Early_active",
                                                    "PM_naive",
                                                    "14dpsham",
                                                    "10dpci_veh",
                                                    "14dpci"))

DimPlot(obj, group.by = "orig.ident")

# Check that it worked as intended
obj@meta.data[["merged_clusters"]] %>% unique()

# Assign identities to cells based on the present cluster
Idents(obj) <- "merged_clusters"

# What levels does the clustering have?
clusters <- levels(obj)

# How many clusters do we have?
cli <- clusters %>% length()
cli

########################## Save the merged cluster object ###################
saveRDS(obj, merged_object)

######################### Check the merged cluster object
obj <- readRDS(merged_object)
DimPlot(obj, label = TRUE)
