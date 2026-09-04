############################### README ###############################
## Author: Inika Prasad, inika.prasad@gmail.com
# Date created: 2025-11-20

# Aims: To obtain classifications for spot identities
# by integrating with single-nuclei transcriptomics dataset
# by label transfer. Colloquially called deconvolution.

# This performs best if you have as similar representation of cell types
# in both the query and reference
# i.e. don't try to map a subset of clusters. It will not do so well.

################################ Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(writexl)
library(readxl)
library(patchwork)

library(ggpointdensity)
library(ggnewscale)

set.seed(100)

source("workflow/functions.R")

#################################### Load variables
whichcondition <- "all"

plot1 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_query_reference.pdf")
plot2 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_deconvolution_UMAPs.pdf")
# plot3 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_predicted_spot_IDs.png")
plot4 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_prediction_scores_violins.pdf")
plot6 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_prediction_scores_ycoord.pdf")
plotx <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_predicted_gene_expressions.pdf")
plot7 <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_", lowerlimit, "_spatial_plots.pdf")
pngfolderhome <- paste0("../results/09_deconvolution/", whichcondition)


output_rds <- paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_seurat_object.rds")

#################################### Load data

# Single-nuclei data
single_nuclei <- readRDS("../../snRNAseq_2025/code/snakemake/results/testerdump/10dpci_ctrl_DRscrub_consolidated.rds")
sn_res <- "new_ident"
sn_reduction <- "umap.harmony"

# Spatial data
spatial <- readRDS("../results/06_dimred_clustering/allsections/allsections_dimred.rds")
spatial_res <- "cca_resolution_1.6"
spatial_reduction <- "integrated.cca"

#################################### subset vehicle only
# whichcondition <- "10dpci_Vehicle"
# single_nuclei <- subset(single_nuclei, subset = orig.ident == "10dpci_veh")
# spatial <- subset(spatial, subset = condition == "vehicle")

# single_nuclei$orig.ident %>% unique()
# spatial$condition %>% unique()

##################################### Initial dataset
# High-contrast color palette for spatial data
spatial_clusters <- spatial@meta.data[[spatial_res]] %>% unique()
selected_cols <- pals::glasbey(30)
names(selected_cols) <- spatial_clusters

pdf(plot1, height = 8, width = 8)
DimPlot(single_nuclei, group.by = sn_res,
        reduction = sn_reduction,
        label = TRUE,
        repel = TRUE,
        pt.size = 0.3)

DimPlot(spatial,
        group.by = spatial_res,
        cols = selected_cols,
        label = TRUE,
        pt.size = 1.5)
dev.off()

#################################### Basic QC of datasets
intersect(rownames(single_nuclei), rownames(spatial)) %>% length()

dim(single_nuclei)
dim(spatial@assays$Spatial)

single_nuclei <- JoinLayers(single_nuclei)

#################################### Label transfer
anchors <- FindTransferAnchors(reference = single_nuclei,
                               query = spatial,
                               reduction = "cca",
                               normalization.method = "LogNormalize")

predictions.assay <- TransferData(anchorset = anchors,
                                  refdata = single_nuclei@meta.data[[sn_res]],
                                  # refdata =  GetAssayData(single_nuclei@assays[['RNA']]),
                                  prediction.assay = TRUE,
                                  weight.reduction = spatial[["integrated.cca"]],
                                  dims = 1:30)

# Assign predictions
spatial[["predictions"]] <- predictions.assay
DefaultAssay(spatial) <- "predictions"

# saveRDS(spatial, file = output_rds)


#################################### Evaluate unmapped clusters
# NOTE:ELEPHANT. START HERE IF YOU'RE NOT RERUNNING THE ANALYSIS
spatial <- readRDS(output_rds)

# If some of the cell types are not found their rowSums will be 0
spatial %>% rowSums() != 0

spatial@assays$predictions$data %>% rowSums()

#################################### Set cutoff for low values
lowerlimit <- 0.01

# Extract assay of prediction values
assaydata <- spatial@assays$predictions$data[, ]

# Initialize datasets where values below the lowerlimit become...
assaydata_cutoff <- assaydata # ...become 0 (for VlnPlots)
assaydataNA <- assaydata # ...become NA (for Spatial & Feature plots)

# Values below lowerlimit become 0
assaydata_cutoff[assaydata_cutoff < lowerlimit] <- 0

# Values below lowerlimit become NA
assaydataNA[assaydataNA < lowerlimit] <- NA

# Quality control for clusters of interest
# NAdf <- assaydataNA["10", ] %>% as.data.frame() %>% rownames_to_column("cell")
# df <- assaydata_cutoff["10", ] %>% as.data.frame() %>% rownames_to_column("cell")
# mappeddf <- full_join(NAdf, df, by = join_by(cell), suffix = c("NAdf", "df"))
# mappeddf %>% View()
# mappeddf %>% filter(is.na(.NAdf) & .df >= lowerlimit) # Should give no results

################################# Set up for plotting (NA values)
# Which clusters do we want plots for
predicted_classes <- sort(rownames(spatial))

# Add the modified data to the spatial_copy object
spatial@assays$predictions$data[, ] <- assaydataNA

################################# Feature Plots of scores
# Visualize cell type distributions on the UMAP
pdf(plot2, width = 25, height = 50)

minvalue <- spatial@assays$predictions$data[, ] %>% min(na.rm = TRUE)

if (minvalue < lowerlimit) {
  stop("There are values below the lower limit,
        ensure data cleaning is done properly.")
} else {
  FeaturePlot(spatial,
              features = predicted_classes,
              keep.scale = "all")  &
    scale_color_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      breaks = c(lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(lowerlimit, 0.25, 0.5, 0.75, 1),
      na.value = "lightgrey"
    )

}

dev.off()


################################## Get predicted IDs
# This plot and proportions are obsolete until I recompute the predicted ID per spot
# based on the mapping of various different clusters

# png(plot3, width = 800, height = 800)
# # Plot predicted IDs for the spatial
# DimPlot(spatial, group.by = "predicted.id", pt.size = 1, label = TRUE)
# dev.off()

# # How many spots in each predicted ID?
# clust_prop <- table(spatial$predicted.id, spatial$condition)
# prop.table(clust_prop, margin = 2) * 100

#################################### Make spatial plots (pdf)
pdf(plot7, height = 20, width = 20)
predicted_classes
# x <- SpatialPlot(object = spatial, features = predicted_classes, ncol = 4, pt.size.factor = 3)

for (i in predicted_classes) {

  minvalue <- spatial@assays$predictions$data[, ] %>% min(na.rm = TRUE)

  if (minvalue < lowerlimit) {
    stop("There are values below the lower limit,
         ensure data cleaning is done properly.")
  } else {

    # Make the spatial plot
    x <- SpatialPlot(object = spatial,
                     features = i,
                     image.alpha = 0.5,
                     ncol = 4,
                     pt.size.factor = 3,
                     stroke = 0.5) &
      scale_fill_gradientn(colors = c("lightgrey",
                                      pals::kovesi.rainbow(256)),
                           values = c(lowerlimit, 1),
                           limits = c(lowerlimit, 1),
                           breaks = c(lowerlimit, 0.5, 1),
                           labels = c(lowerlimit, 0.5, 1)) &
    theme(
      legend.spacing = unit(1, "cm")
    )
    print(x)

  }



}

dev.off()

#################################### Make spatial plots (pngs for publication)
# Specify the section for which you want plots
image <- "vehicle.5"
image <- "vehicle.7"
image <- "c.CpE.6"
image <- "naive.1"

# Print a plot for every predicted class
for (i in predicted_classes) {

  minvalue <- spatial@assays$predictions$data[, ] %>% min(na.rm = TRUE)

  if (minvalue < lowerlimit) {
    stop("There are values below the lower limit,
         ensure data cleaning is done properly.")
  }

  # Make a filepath for the image
  filename <- file.path(pngfolderhome,
                        image,
                        paste0(image, "_", i, "_NA", lowerlimit, ".png"))
  filename

  # Make a filepath for the image
  filename_legend <- file.path(pngfolderhome,
                               image,
                               paste0(image, "_", i, "_NA", lowerlimit, "_legend", ".png"))
  filename_legend

  # Create the directory to store images if it doesn't exist already
  if (!dir.exists(dirname(filename))) {
    dir.create(dirname(filename))
  }

  x <- SpatialPlot(
    object = spatial,
    images = image,
    features = i,
    image.alpha = 0.5,
    pt.size.factor = 3.5,
    stroke = 0.5
  ) &
    scale_fill_gradientn(
      colors = pals::kovesi.rainbow(256),
      limits = c(lowerlimit, 1),
      na.value = "lightgrey",
      breaks = c(0, lowerlimit, 0.25, 0.5, 0.75, 1),
      labels = c(0, lowerlimit, 0.25, 0.5, 0.75, 1)
    ) &
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      # make legend bigger
      legend.key.height = unit(3, "cm"),
      legend.key.width  = unit(0.6, "cm"),
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 30)
    )

  # Extract plot and legend separately
  plot_no_legend <- x +
    theme(legend.position = "none") +
    ggtitle(paste0(i, ", NA threshold: ", lowerlimit))
  legend <- cowplot::get_legend(x)

  # Save the plot and legend
  ggsave(filename, plot_no_legend)
  ggsave(filename_legend, legend)


}

################################# Set up for plotting (0 values)
# Add the modified data to the spatial_copy object
spatial@assays$predictions$data[, ] <- assaydataNA

################################## Get Violin plot of values
sum(is.na(spatial[[predicted_classes[1]]]))

# Making NA values causes problems for Violin Plot
spatial@assays$predictions$data[, ] <- assaydatacutoff

pdf(plot4, width = 25, height = 35)
VlnPlot(spatial,
        features = predicted_classes,
        alpha = 0.02,
        ncol = 4,
        cols = selected_cols,
        group.by = spatial_res)
dev.off()

####################### Quantile plots
pdf(paste0("../results/09_deconvolution/", whichcondition, "/", whichcondition, "_percentiles.pdf"), height = 2, width = 3)

# SHOULD THIS BE DONE AFTER CHANGING 0s to NA OR BEFORE??
probs  <- seq(0, 1, by = 0.01) # look at every 1percentile
qvals  <- quantile(assaydata, probs = probs, na.rm = TRUE)

df <- data.frame(
  quantile = probs,
  value = as.numeric(qvals)
)

quantplot <- ggplot(df, aes(x = quantile, y = value)) +
  geom_line() +
  # scale_y_log10() +
  labs(x = "Quantile", y = "Value", title = i) +
  theme_minimal()

print(quantplot)

dev.off()

#################################### Evaluate predicted expression of genes

predictions.assay.genes <- TransferData(anchorset = anchors,
                                  # refdata = single_nuclei@meta.data[[sn_res]],
                                  refdata =  GetAssayData(single_nuclei@assays[['RNA']]),
                                  prediction.assay = TRUE,
                                  weight.reduction = spatial[["integrated.cca"]],
                                  dims = 1:30)

# Assign predictions
spatial[["predictions.genes"]] <- predictions.assay.genes
DefaultAssay(spatial) <- "predictions.genes"

genelist <- read_xlsx("../results/09_deconvolution/epipeeloutgenes.xlsx") %>% pull(.)

pdf(plotx, height = 12)
# Visualize gene expression: predicted and true
for (gene in genelist) {

  DefaultAssay(spatial) <- "Spatial"
  if (gene %in% rownames(spatial)) {
    p1 <- FeaturePlot(spatial, features = gene, pt.size = 0.6, order = TRUE) +
    ggtitle(paste0("true_", gene))
  } else {p1 <- FeaturePlot(spatial, "nFeature_Spatial", pt.size = 0.6, cols = c("grey", "green"))}

  DefaultAssay(spatial) <- "predictions.genes"
  if (gene %in% rownames(spatial)) {
  p2 <- FeaturePlot(spatial, features = gene, pt.size = 0.6, order = TRUE) +
    ggtitle(paste0("predicted_", gene))
  } else {p2 <- FeaturePlot(spatial, "nFeature_Spatial", pt.size = 0.6, cols = c("grey", "green"))}

  print(p1 + p2)

}
dev.off()

########################## Visualize localization & scoring: base to apex
# Extract coordinate information from metadata
meta <- spatial@meta.data

# Group by sample and get normalized x and y coordinates
meta <- meta %>%
  group_by(sample) %>%
  mutate(
    x_norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
    y_norm = (y - min(y, na.rm = TRUE)) / (max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  ) %>%
  ungroup()


pdf(plot6)
# Choose a cell type of interest
pred_scores <- colnames(meta)[grep("prediction", colnames(meta))] %>% sort()

for (pred_score in pred_scores) {

  myscore <- pred_score
  print(myscore)

  # Remove spots with no scores
  meta2 <- meta %>% filter(!!sym(myscore) != 0)

  # Pretty density plot
  p <- ggplot(meta2, aes(x = !!sym(myscore), y = y_norm)) +
    geom_density_2d_filled() +
    geom_point(alpha = 0.3) +
    theme_minimal() +
    theme(text = element_text(size = 20))

  print(p)

}

dev.off()