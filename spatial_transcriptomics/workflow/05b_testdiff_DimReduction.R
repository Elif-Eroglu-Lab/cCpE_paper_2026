# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: For a Seurat object, test a range of different dimensions for RunUMAP
# This is useful for determining which dimensions to proceed with.

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

###############################  Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(clustree)
library(gridExtra)
library(grid)
library(patchwork)

set.seed(100)

########################### Load tester variables
# rds <- obj
# ndims <- 15
# reduction2test <- "integrated.cca"
# # txt1 <- 
# plot1 <- paste0("diff_ndims_", reduction2test, ".pdf")
# plot4 <- snakemake@output[["plot4"]]

########################### Load snakemake variables
rds <- snakemake@input[["input_rds"]]
ndims <- as.numeric(snakemake@params$ndims) # User-defined ndims to test
reduction2test <- snakemake@params$reduction2test # User-defined dimensionality reduction

txt1 <- snakemake@output[["txt1"]]
plot1 <- snakemake@output[["plot1"]]
plot4 <- snakemake@output[["plot4"]]

########################### Load data
sfi <- readRDS(rds)

########################### Test Snakemake variables
# Ensure the requested reduction exists in the object
if (!reduction2test %in% names(sfi@reductions)) {
  stop(paste0(
    reduction2test,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(sfi@reductions), collapse = ", ")
  ))
}

########################### Specify the dimensions to test
diff_ndims <- c(ndims-10, ndims-5, ndims, ndims+5, ndims+10, ndims+15)

diff_ndims <- diff_ndims[diff_ndims >= 0]
diff_dimred_names <- paste0("UMAP_", reduction2test, diff_ndims)

########################### Record testing info
# Record the ndims and reduction used in an accompanying text file
writeLines(c(paste0("reduction used for UMAP:", reduction2test),
             paste0("ndims used: ", ndims),
             "Resulting reductions:",
             diff_dimred_names),
           con = txt1)

sfi@meta.data %>% colnames()
sfi@active.ident %>% levels()

########################### Test different dimensionality reductions

Idents(sfi) <- sfi$orig.ident # So that the UMAPs are colored by orig.ident

# UMAP dimensionality reduction & basic plotting
pdf(plot1, width = 20)
for (i in seq_along(diff_ndims)) {
  ndim2try <- diff_ndims[i]
  dimred_name <- diff_dimred_names[i]
  sfi <- RunUMAP(sfi,
    dims = 1:ndim2try,
    reduction = reduction2test, # name of input reduction, defined by user
    reduction.name = dimred_name # name of output UMAP
  )

  p1 <- DimPlot(sfi,
                reduction = dimred_name) +
    ggtitle(dimred_name)

  p2 <- DimPlot(sfi,
                reduction = dimred_name,
                group.by = "condition") +
    ggtitle(dimred_name)

  p3 <- DimPlot(sfi,
                reduction = dimred_name,
                group.by = "batch",
                cols = c("purple", "darkgreen")) +
    ggtitle(dimred_name)

  print(p1 + p2 + p3)
}
dev.off()

##########################

# Plot nFeature & nCount for each UMAP
pdf(plot4, height = 14)
for (i in seq_along(diff_ndims)) {
  dimred_name <- diff_dimred_names[i]
  p4 <- FeaturePlot(sfi,
                    reduction = dimred_name,
                    features = c("nCount_Spatial", "nFeature_Spatial"),
                    label.size = 2)
  # p5 <- FeaturePlot(sfi,
  #                   reduction = dimred_name,
  #                   features = c("filt_nCount_Spatial", "filt_nFeature_Spatial"),
  #                   label.size = 2)
  print(p4)
  # print(p5)
}
dev.off()

###############################  Close log file
sink(type = "message")
sink()
close(logfile)