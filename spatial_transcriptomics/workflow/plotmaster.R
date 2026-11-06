############################### README ###############################
# Aim: Takes an excel sheet as input and creates a custom plot grid
# Input: Excel sheet with columns group.by, group, features
# Output: custom_plots.pdf where each row is a section,
# and each column is a grouping (like epicardial spots) or a feature (like CLD6)
# Author: Inika Prasad, inika.prasad@ki.se
# Created: 05 September, 2025

###############################  Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(readxl)
library(pals) #conda install r-pals
library(viridis)

set.seed(100)

source("workflow/functions.R")

############################ Setup object ###########################
# Load data
obj <- readRDS("../results/06_dimred_clustering/allsections/allsections_dimred.rds")
myclusters <- "cca_resolution_1.6"
Idents(obj) <- myclusters

# High-contrast color palette
# selected_cols <- pals::glasbey(30)
# names(selected_cols) <- Idents(obj) %>% levels()
# show_cols(selected_cols)

# Manually chosen color palette
selected_cols <- c(
"01" = "#EC8398",
"02" = "#D62246", #replaced
"03" = "#F4B8C4",
"04" = "#7C1328",
"05" = "#FAA916",
"06" = "#48B8D0",
"07" = "#267E92",
"08" = "#CE7FF0",
"09" = "#FF495C",
"10" = "#00FFBE", #BZ done
"11" = "#62466B",
"12" = "#D741A7",
"13" = "#D89D6A",
"14" = "#6C8CB2",
"15" = "#B1CC71",
"16" = "#FE8F42",
"17" = "#799636",
"18" = "#FFD300",
"19" = "#720055",
"20" = "#b23fd8",
"21" = "#02AD24",
"22" = "#C8FF00")

levels(obj$animal_readable)
obj$animal_readable <- factor(obj$animal_readable,
                              levels = c("naive_A",
                                         "vehicle_A", "vehicle_B", "vehicle_C",
                                         "c-CpE_A", "c-CpE_B", "c-CpE_C"))

############################ Make the legend #######################
# Convert to data frame
variable <- "batch"
variable <- "animal_readable"
variable <- "sample"
nlevels <- length(levels(obj@meta.data[[variable]]))
selected_cols <- pal_hue()(nlevels)
names(selected_cols) <- levels(obj@meta.data[[variable]])

df <- data.frame(
  name = names(selected_cols),
  color = unname(selected_cols),
  x = seq_along(selected_cols),
  y = 1
)

ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_tile(width = 0.9, height = 0.3) +  # skinny rectangles
  scale_fill_identity() +                 # use actual color values
  # geom_text(aes(label = name), vjust = -1, size = 3) +
  theme_void() +
  theme(legend.position = "none")

filename <- paste0("../results/figures_spatial/ST_clusters_legend_", variable, ".png")
ggsave(filename, height = 1, width = nlevels)
ggsave("../results/figures_spatial/ST_clusters_legend_withnames.png", height = 1, width = 25)

############################# Plot basic UMAPs #####################
ptsize <- 1.2

DimPlot(obj, group.by = myclusters, cols = selected_cols, label = FALSE, pt.size = ptsize)
ggsave("../results/figures_spatial/ST_dimplot.png")

DimPlot(obj, group.by = "sample", pt.size = ptsize)
ggsave("../results/figures_spatial/ST_dimplot_samples.png")

DimPlot(obj, group.by = "animal_readable", pt.size = ptsize)
ggsave("../results/figures_spatial/ST_dimplot_animal.png")

DimPlot(obj, group.by = "batch", pt.size = ptsize)
ggsave("../results/figures_spatial/ST_dimplot_batch.png")

DimPlot(obj, group.by = myclusters, split.by = "condition",  cols = selected_cols, pt.size = ptsize, ncol = 1)
ggsave("../results/figures_spatial/ST_dimplot_conditionsplit.png", height = 25)

############################ Plot QC UMAPs ########################
joined_obj <- JoinLayers(obj)
# Compute new values of nFeature and nCount after filtering
joined_obj$filt_nCount_Spatial <- colSums(joined_obj@assays$Spatial$counts)
joined_obj$filt_nFeature_Spatial <- colSums(joined_obj@assays$Spatial$counts != 0)

FeaturePlot(joined_obj, features = "nFeature_Spatial", pt.size = 1.2) +
FeaturePlot(joined_obj, features = "filt_nFeature_Spatial", pt.size = 1.2) &
    scale_colour_viridis(option = "viridis")

FeaturePlot(joined_obj, features = "nCount_Spatial", pt.size = 1.2) +
FeaturePlot(joined_obj, features = "filt_nCount_Spatial", pt.size = 1.2) &
    scale_colour_viridis(option = "viridis")

# it's functionally equivalent
FeaturePlot(joined_obj, features = "nCount_Spatial", pt.size = 1.2) &
    scale_colour_viridis(option = "viridis")
ggsave("../results/figures_spatial/ST_nCount.png")
ggsave("../results/figures_spatial/ST_nCount.pdf")


FeaturePlot(joined_obj, features = "nFeature_Spatial", pt.size = 1.2) &
    scale_colour_viridis(option = "viridis")
ggsave("../results/figures_spatial/ST_nFeature.png")
ggsave("../results/figures_spatial/ST_nFeature.pdf")


# VlnPlot(joined_obj, features = "nCount_Spatial", alpha = 0.5, group.by = "animal_readable")
VlnPlot(joined_obj, features = "nCount_Spatial", alpha = 0, group.by = "sample")
ggsave("../results/figures_spatial/ST_nCount_VlnPlot_nopoints.pdf", height = 7)

# VlnPlot(joined_obj, features = "nCount_Spatial", alpha = 0.5, group.by = "animal_readable")
VlnPlot(joined_obj, features = "nFeature_Spatial", alpha = 0.1, group.by = "sample")
ggsave("../results/figures_spatial/ST_nFeature_VlnPlot.pdf", height = 7)

############################# Plot Border Zone #####################
# Set everything as grey except cluster(s) of interest
selected_cols <- rep("darkgrey", 30)
names(selected_cols) <- Idents(obj) %>% levels()
selected_cols["10"] <- "#00FFBE"


x <- plotall(obj,
  crop = TRUE,
  # selected_clustering = myclusters,
  selected_cols = selected_cols,
  image.alpha = 0,
  pt.size.factor = 3.2,
  # cells.highlight = cells
)
x$naive
x$vehicle.5
x$c.CpE.6
ggsave("../results/figures_spatial/CHANGETHIS.png")

# Save all the images
pdf("../results/figures_spatial/BZonly.pdf")
x
dev.off()

############################## Make dotplot
{ # Top markers for the PRE-MERGING Epi peelout
markers <- read_xlsx("../results/07_clustermarkers/allsections/allsections_cca_resolution_1.6_markers.xlsx")
markers$cluster <- factor(markers$cluster, levels = names(selected_cols))
markers$score <- as.numeric(markers$score)

# Map back to original identifiers
identifiers <- read_xlsx("../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.xlsx")
identifiermap <- setNames(identifiers$gene_id, identifiers$genename_Hs_PwNameFormat)

# Make supplementary file for Spatial markers
colnames(markers)
markers_supp <- markers %>%
  dplyr::select(-score2) %>%
  dplyr::mutate(gene_id = identifiermap[gene]) %>%
  dplyr::relocate(gene_id, gene,cluster)
write_xlsx(markers_supp, "../results/figures_spatial/SpatialClusters_markers.xlsx")


View(markers)

# Select markers to plot
genes2plot <- markers %>%
  group_by(cluster) %>%
  arrange(cluster, -score, , -(abs(avg_log2FC))) %>%
  slice_head(n = 4) %>%
  pull(gene) %>%
  unique()
  
  # %>%
  # rev() # for vertical plot

width <- 0.15 * length(genes2plot) + 7
height <- 6

label <- "ST_4markers_byscore_and_log2FC_SCALED"
}

p <- DotPlot(obj,
        scale = TRUE,
        unique(genes2plot),
        group.by = myclusters) +
  scale_colour_viridis(option = "viridis") +
  scale_x_discrete(labels = function(x) stringr::str_to_sentence(identifiermap[x])) +
  RotatedAxis()
  # +
  # coord_flip()
p

ggsave(paste0("../results/figures_spatial/HORIZONTAL", label, "dotplot.pdf"),
       width = width,
       height = height,
       units = "in")
