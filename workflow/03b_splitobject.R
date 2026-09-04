############################## README
# Author: Inika Prasad
# Aim: To split the Seurat object by the variable of interest.
# E.g. each layer corresponds to a section (default) or
# each layer corresponds to a batch (correct for batch effect)

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

###############################  Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

source("workflow/functions.R")

set.seed(100)

############################### Tester variables
# input_rds <- "../results/03_Filtering/allsections/allsections_filtered.rds"
# split_variable <- "batch"

############################### Load Snakemake variables
input_rds <- snakemake@input[["input_rds"]]
output_rds <- snakemake@output[["rds"]]

split_variable <- snakemake@params[["split_variable"]]

############################### Load data
obj <- readRDS(input_rds)

############################### Split data as needed
obj <- JoinLayers(obj)

obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj@meta.data[[split_variable]])

obj
print(paste0("Object split into layers based on: ", split_variable))

############################### Save data

saveRDS(obj, file = output_rds)