# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: For a Seurat object, test a range of different dimensions for RunUMAP
# This is useful for determining which dimensions to proceed with.

############################# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################# Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(clustree)
library(gridExtra)
library(grid)
library(patchwork)

set.seed(100)

########################### Load snakemake variables
rds <- snakemake@input[["rds"]]
input_seurat <- snakemake@params$input_seurat
ndims <- as.numeric(snakemake@params$ndims) # User-defined ndims to test
reduction <- snakemake@params$reduction # User-defined dimensionality reduction
knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col
join_layers <- snakemake@params$join_layers # Should be "TRUE" or "FALSE"
plot_doublets <- snakemake@params$plot_doublets # Should be "TRUE" or "FALSE"

txt1 <- snakemake@output[["txt1"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]

# output_rds <- snakemake@output[["output_rds"]]

print(knownmarker_file)
print(as.character(input_seurat))

if(is.null(input_seurat)){
  stop(paste0("sample name is missing, and is needed for dotplot labeling"))
}
########################### Load data
sfi <- readRDS(rds)

## If we are working with integrated data, then we need to JoinLayers
## So that the feature plots don't plot NA values as dark grey
## But instead plot them as light grey values of
if (join_layers == "TRUE") {
  sfi <- JoinLayers(sfi)
}

# Load list of marker genes
goi <- read.csv(knownmarker_file, sep = ";")

########################### Test Snakemake variables
# Ensure the requested reduction exists in the object
if (!reduction %in% names(sfi@reductions)) {
  stop(paste0(
    reduction,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(sfi@reductions), collapse = ", ")
  ))
}

# Ensure marker file exists
if (!file.exists(knownmarker_file)) {
  stop(paste0(knownmarker_file, "does not exist. Please check path"))
}

# Ensure marker file has the necessary column
if (!knownmarker_id_col %in% colnames(goi)) {
  stop(paste0(
    knownmarker_id_col,
    " is not a column in the marker file.\n",
    "Available columns are: ",
    paste(colnames(goi), collapse = ", ")
  ))
}

############################## Place samples in an intuitive order
sfi$orig.ident <- factor(sfi$orig.ident,
                         levels = c("Early_active",
                                    "PM_naive",
                                    "14dpsham",
                                    "10dpci_veh",
                                    "14dpci",
                                    "10dpci_cCpE"))

Idents(sfi) <- sfi$orig.ident

########################### Define dimensions to test

diff_ndims <- c(ndims-10, ndims-5, ndims, ndims+5, ndims+10, ndims+15)

diff_ndims <- diff_ndims[diff_ndims >= 0]
diff_dimred_names <- paste0("UMAP_", reduction, diff_ndims)

# Record the ndims and reduction used in an accompanying text file
writeLines(c(paste0("reduction used for UMAP:", reduction),
             paste0("ndims used: ", ndims),
             "Resulting reductions:",
             diff_dimred_names),
           con = txt1)

sfi@meta.data %>% colnames()
sfi@active.ident %>% levels()
Idents(sfi) <- sfi$orig.ident # So that the UMAPs are colored by orig.ident

# UMAP dimensionality reduction & basic plotting
pdf(plot1)
for (i in seq_along(diff_ndims)) {
  ndim2try <- diff_ndims[i]
  dimred_name <- diff_dimred_names[i]
  sfi <- RunUMAP(sfi,
    dims = 1:ndim2try,
    reduction = reduction, # name of input reduction, defined by user
    reduction.name = dimred_name # name of output UMAP
  )

  p <- DimPlot(sfi,
               reduction = dimred_name) +
    ggtitle(dimred_name)

  print(p)
}
dev.off()

if (plot_doublets == "TRUE") {

  # Plot doublet confidence for each UMAP
  pdf(plot3, width = 14)
  for (i in seq_along(diff_ndims)) {
    dimred_name <- diff_dimred_names[i]
    p3 <- DimPlot(sfi,
                  reduction = dimred_name,
                  group.by = "scrublet_DropletType") +
      theme(legend.position = "bottom")

    p3_2 <- FeaturePlot(sfi,
                        reduction = dimred_name,
                        features = c("scrublet_Scores"),
                        cols = c("grey", "darkgreen"),
                        label.size = 2)
    print(p3 + p3_2)
  }
  dev.off()

  # Plot nFeature & nCount for each UMAP
  pdf(plot4, width = 14)
  for (i in seq_along(diff_ndims)) {
    dimred_name <- diff_dimred_names[i]
    p4 <- FeaturePlot(sfi,
                      reduction = dimred_name,
                      features = c("nCount_RNA", "nFeature_RNA"),
                      label.size = 2)

    print(p4)
  }
  dev.off()

}

# Evaluate Cell Cycle
pdf(plot5, width = 21)
for (i in seq_along(diff_ndims)) {
  dimred_name <- diff_dimred_names[i]

  p5_1 <- DimPlot(sfi,
                  reduction = dimred_name,
                  group.by = "Phase") +
  theme(legend.position = "bottom")

  p5_2 <- FeaturePlot(sfi,
                      reduction = dimred_name,
                      features = "S.Score",
                      cols = c("grey", "#619CFF"),
                      label.size = 2)

  p5_3 <- FeaturePlot(sfi,
                      reduction = dimred_name,
                      features = "G2M.Score",
                      cols = c("grey", "#00BA38"),
                      label.size = 2)

  print(p5_1 + p5_2 + p5_3)
}
dev.off()

# Evaluate sample split
pdf(plot6, width = 7 * length(unique(sfi$orig.ident)))
for (i in seq_along(diff_ndims)) {
  dimred_name <- diff_dimred_names[i]
  p6 <- DimPlot(sfi,
                reduction = dimred_name,
                split.by = "orig.ident")
  print(p6)
}
dev.off()

########################################### Plot basic feature plots
# Specify genes of interest (known marker genes)
print(goi)
# Clean up list
goi <- goi[!is.na(goi[[knownmarker_id_col]]) & goi[[knownmarker_id_col]] != "", ]
goi <- goi[!duplicated(goi[[knownmarker_id_col]]), ]

names_goi <- goi[[knownmarker_id_col]]
goi2plot <- names_goi[names_goi %in% rownames(sfi)]
if (length(goi2plot) == 0) {
  warning("No marker genes from the marker file are present in the Seurat object. Please check your marker file and Seurat object rownames.")
  stop("No marker genes found in Seurat object. Exiting.")
}

chunk_size <- 6  # max plots per page
gene_chunks <- split(goi2plot, ceiling(seq_along(goi2plot) / chunk_size))

# UMAP dimensionality reduction
pdf(plot2, width = 11.7, height = 8.3)
for (i in seq_along(diff_ndims)) {
  dimred_name <- diff_dimred_names[i]

  for (x in seq_along(gene_chunks)) {
    genes <- gene_chunks[[x]]

    # One FeaturePlot per gene
    plots <- lapply(genes, function(g) {
      FeaturePlot(
        sfi,
        features = g,
        reduction = dimred_name,
        raster = TRUE
      )
    })

    # Arrange in grid with a title
    p <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = paste0("Feature Plots (Page", x, "): ", input_seurat)
      )

    print(p)
  }
}
dev.off()

############################# # Save the dimensionality reduced data
# saveRDS(sfi, output_rds)

############################# Close log file
sink(type = "message")
sink()
close(logfile)