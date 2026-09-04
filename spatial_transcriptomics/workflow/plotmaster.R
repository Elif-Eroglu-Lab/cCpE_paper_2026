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
markers$score2 <- as.numeric(markers$score2)

View(markers)

# Select markers to plot
genes2plot <- markers %>%
  group_by(cluster) %>%
  arrange(cluster, -score, , -(abs(avg_log2FC))) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()
  
  # %>%
  # rev() # for vertical plot

label <- "ST_5markers_byscore_and_log2FC_SCALED"
}

p <- DotPlot(obj,
        scale = TRUE,
        unique(genes2plot),
        group.by = myclusters) +
  scale_colour_viridis(option = "viridis") +
  scale_x_discrete(labels = function(x) stringr::str_to_sentence(x)) +
  RotatedAxis()
  # +
  # coord_flip()
p

ggsave(paste0("../results/figures_spatial/HORIZONTAL", label, "dotplot.pdf"),
      #  height = 0.15 * length(genes2plot) + 1.8,
      #  width = 10,
       width = 0.15 * length(genes2plot) + 7,
       height = 8,
       units = "in")

############################ Setup stuff to plot ###########################

# Make sure that the groups actually exist in the group.by dataset!!
stuff2plot <- read_xlsx("misc_totrash/custom_plot_grids.xlsx") %>% as.data.frame()
stuff2plot[is.na(stuff2plot)] <- ""
stuff2plot

DimPlot(sfi, group.by = clustering, label = TRUE)


FeaturePlot(obj, toplot, pt.size = 1.2)
SpatialFeaturePlot(sfi, "BNC1")

toplot <- c("BNC1", "ATP6V1B1", "SLC43A2-LOC138304169", "ATP6V1C2", "NEURL1", "GRPR", "ATP6V0A4-LOC138287520", "SGK2", "TBC1D24-LOC138246090", "LOC138300737", "ATP6AP1-LOC138288350", "TACR1", "PGR", "ATP6V1A", "TMEM213", "TFCP2L1-PW", "GNAZ", "RBM47", "HEPACAM2-PW", "RAPGEF3", "SPTBN2", "TFCP2L1-LOC138299295", "PDE8B", "UGT3A1-LOC138282918", "SLC10A4", "PLXDC2", "GALR1-PW", "BTN1A1-LOC138299205", "PCSK6", "LOC138246465", "TMEM171", "ATP6V0D2", "ATP6V1G3", "CFTR", "FYB2", "MUC1-PW", "CLCNKB-LOC138300562", "HSD11B2", "PRSS33-LOC138304300", "LPAR3-PW", "BAIAP2L1", "GJB7", "LOC138282493")
toplot <- "BNC1"
toplot <- c("LOC138301428", "LOC138267955", "LOC138260960", "LOC138294116")

x <- plotall(obj,
        crop = TRUE,
        features2plot = toplot,
        image.alpha = 0,
        pt.size.factor = 3)

pdf("LOCgenes.pdf", width = 18, height = 12)
x
dev.off()

############################ Setup stuff to plot ###########################

# Turn empty values to NULL for the plotmaster function
df2 <- stuff2plot %>%
  mutate(across(everything(), ~ lapply(.x, function(x) if (x == "") NULL else x)))

df2


# Modified function
# Function: Plot clusters in each image
plotmaster <- function(object,
                       features = NULL,
                       group.by = NULL,
                       cols = NULL,
                       image.alpha = 0,
                       legend.dot.size = 7,
                       label = FALSE,
                       ...) {

  images <- names(object@images)

  # Decide what we’re plotting
  if (!is.null(features) && !is.null(group.by)) {
    stop("Please provide either 'features' or 'group.by', not both.")
  } else if (!is.null(features)) {
    mode <- "Feature"
  } else if (!is.null(group.by)) {
    mode <- "Group"
  } else {
    stop("You must provide either 'features' or 'group.by'.")
  }

  plots <- lapply(images, function(s) {
    p <- SpatialPlot(object = object,
                     features = features,
                     group.by = group.by,
                     cols = cols,
                     images = s, # Specify the image name
                     image.alpha = image.alpha,
                     alpha = c(0.9, 1), # Transparency for image and spots
                     label = label,
                     ...)
      # Remove title of plots for this application
      # ggtitle(paste(d$sample[d$id == s]))
      # ggtitle(s)

    if (mode == "Feature") {
      # Find maximum value for scale
      df <- FetchData(object = object, vars = features)
      maxval <- round(x = max(df))

      p <- p + scale_fill_gradientn(limits = c(0, maxval),
                                    # Grey for 0, Viridis for >0
                                    colors = c("grey", viridis::plasma(100))) +
        theme(legend.position = "right",
              legend.key.size = unit(0.5, "cm"),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9))

    } else if (mode == "Group") {
    #   p <- p +
    #     guides(fill = "none")
      p <- p + ggtitle(group.by)
    }
    return(p)
  })

  names(plots) <- images

  return(plots)
}


########################## Make plots ###############################
obj <- sfi

# Initialise list of plots
p <- list()

# For every entry, make the plots
for (i in (1:nrow(df2))) {

  # Identify which group to highlight (may be NULL)
  highlight <- df2$group[[i]]

  # Get your spots highlighted if there is a group to highlight
  if (!is.null(highlight)) {
    selected_cols <- selected_cols_blank
    selected_cols[[highlight]] <- bgcols[[highlight]]
  } else {
    selected_cols <- mypalette
    names(selected_cols) <- obj@meta.data[[df2$group.by[[i]]]] %>% levels()
  }

  # Make plots & add them to a list
  p[[i]] <- plotmaster(object = obj,
                       image.alpha = 0.3,
                       pt.size.factor = 2,
                       group.by = df2$group.by[[i]],
                       features = df2$features[[i]],
                       cols = selected_cols)

}
######################## Combine & Save plots #######################
plotlist <- unlist(p, recursive = FALSE)

plotlist[[1]]
length(plotlist)

# # Make
# pdf("vehicle_8_cluster_markers_spatial.pdf",
#     height = 12,
#     width = 25,
#     bg = "black")
# dev.off()

pdf("clusterings2.pdf", width = 28, height = 35)
patchwork::wrap_plots(plotlist,
            ncol = nrow(df2),
            nrow = length(obj@images),
            byrow = FALSE) & theme(plot.margin = margin(0.02, 0.02, 0.02, 0.02))
dev.off()

##########

selected_cols <- mypalette
names(selected_cols) <- obj@meta.data[["cca_resolution_1"]] %>% levels()

selected_cols

# Baseline grey color for all clusters
selected_cols_blank <- rep("grey", 20)
names(selected_cols_blank) <- obj@meta.data[["cca_resolution_1"]] %>% levels()
# selected_cols_blank[["13"]] <- "magenta"
selected_cols_blank[["01"]] <- "blue"
selected_cols_blank[["02"]] <- "purple"

selected_cols_blank

x <- plotall(obj,
             selected_clustering = "cca_resolution_1",
             selected_cols = selected_cols_blank,
             crop = TRUE,
             image.alpha = 0,
             pt.size.factor = 2)

pdf("clust_res1_musclestomerge.pdf")
x
dev.off()


BM_vsremote <- read_xlsx("../results/08_DGE/ZAVIJAVA_DGE_results.xlsx", sheet = "sig_degs")
BMveh_vsBMccpe <- read_xlsx("../results/08c_refinedDEGS/ZOSMA_ZIBAL_refined.xlsx")

intersect(BM_vsremote$gene, BMveh_vsBMccpe$gene)
