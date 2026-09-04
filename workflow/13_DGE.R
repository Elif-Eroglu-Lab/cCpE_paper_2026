############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

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

################################# Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

unique_code <- snakemake@params$unique_code
clustering <- snakemake@params$clustering
reduction <- snakemake@params$reduction

ds <- snakemake@params$downsampling
covariates <- snakemake@params$covariates
test_use <- snakemake@params$test_use
assay <- snakemake@params$assay

# For individual DGE
tst <- snakemake@params$tst
ref <- snakemake@params$ref
grp <- snakemake@params$grp

# For mass DGE
tst_grp <- snakemake@params$tst_grp
ref_grp <- snakemake@params$ref_grp
grping <- snakemake@params$grping
clusterwise <- snakemake@params$clusterwise

# Outputs
outfile <- snakemake@output[["outfile"]]
outplot <- snakemake@output[["outplot"]]
outres <- snakemake@output[["outres"]]

################################# Basic info

# Helper function to format values cleanly
fmt <- function(x) {
  if (is.null(x)) {
    "NULL"
  } else if (length(x) > 1) {
    paste(x, collapse = ", ")
  } else {
    as.character(x)
  }
}

# Write file
writeLines(
  c(
    paste0('input_rds: "', fmt(input_rds), '"'),
    paste0('unique_code: "', fmt(unique_code), '"'),
    paste0('clustering: "', fmt(clustering), '"'),
    paste0('reduction: "', fmt(reduction), '"'),
    paste0('downsampling: "', fmt(ds), '"'),
    paste0('covariates: "', fmt(covariates), '"'),
    paste0('test_use: "', fmt(test_use), '"'),
    paste0('assay: "', fmt(assay), '"'),

    "",
    "# Individual DGE",
    paste0('tst: "', fmt(tst), '"'),
    paste0('ref: "', fmt(ref), '"'),
    paste0('grp: "', fmt(grp), '"'),

    "",
    "# Mass DGE",
    paste0('tst_grp: "', fmt(tst_grp), '"'),
    paste0('ref_grp: "', fmt(ref_grp), '"'),
    paste0('grping: "', fmt(grping), '"'),
    paste0('clusterwise: "', fmt(clusterwise), '"')
  ),
  con = outfile
)

################################# Load dataset
obj <- readRDS(input_rds)

# Need to ensure this is already done.
obj2 <- JoinLayers(obj)

# Save UMAP for quick recall of which object this is
DimPlot(obj2, group.by = clustering, reduction = reduction, pt.size = 1, label = TRUE)
ggsave(outplot)

# Source the file where the metadata columns are created
source("workflow/08b_make_metadata4DEGs.R")

#  Source the functions
source("workflow/functions.R")

################################ Determine pairs to compare

summary <- list()

# If mass dge is to be run:
if (clusterwise == "TRUE") {

  # Extract metadata
  meta <- obj2@meta.data

  comp_table <- table(meta[[grping]], meta[[clustering]]) %>%
    as.data.frame.matrix() %>%
    # Remove conditions that are not of interest
    filter(rownames(.) %in% c(ref_grp, tst_grp)) %>%
    # Order rownames to match ref_grp, tst_grp
    slice(match(c(ref_grp, tst_grp), rownames(.))) %>%
    # Remove the clusters with too few cells in either condition
    select(where(is.numeric)) %>%
    select(where(~ all(.x >= 5, na.rm = TRUE)))

  print(comp_table)

  # Make pairs for DGE comparisons
  pairs <- lapply(colnames(comp_table),
                  function(cl) {
                                paste0(rownames(comp_table), "_", cl)})

  # Name pairs for unique_code
  pairs <- setNames(pairs, paste0(unique_code, colnames(comp_table)))
  pairs

  for (i in seq_along(pairs)) {

    unique_subcode <- names(pairs[i])
    ref <- pairs[[i]][1]
    tst <- pairs[[i]][2]
    grp <- grp

    # Outputs
    table1 <- paste0("../results/08_DGE/", unique_code, "/", unique_subcode, "_DGE_results.xlsx")
    plot1 <- paste0("../results/08_DGE/", unique_code, "/", unique_subcode, "_DGE_results.pdf")
    plot2 <- paste0("../results/08_DGE/", unique_code, "/", unique_subcode, "_DEG_expression_clusters.pdf")

    # Run the dge_function
    dge_function(seurat_object = obj2,
                 unique_code = unique_subcode,
                 reduction = reduction,
                 tst = tst,
                 ref = ref,
                 grp = grp,
                 ds = ds,
                 covariates = covariates,
                 test_use = test_use,
                 assay = assay,
                 table1 = table1,
                 plot1 = plot1,
                 plot2 = plot2)

    # Extract summary of the run
    summary_x <- read_xlsx(table1, sheet = "run_details")
    summary[[i]] <- summary_x
  }

  # Save summary pf the run
  summary_df <- do.call(rbind, summary)
  write_xlsx(summary_df,
             outres)

} else { # FOR INDIVIDUAL DGE

  # Outputs
  table1 <- paste0("../results/08_DGE/", unique_code, "/", unique_code, "_DGE_results.xlsx")
  plot1 <- paste0("../results/08_DGE/", unique_code, "/", unique_code, "_DGE_results.pdf")
  plot2 <- paste0("../results/08_DGE/", unique_code, "/", unique_code, "_DEG_expression_clusters.pdf")

  # Run the dge_function
  dge_function(seurat_object = obj2,
               unique_code = unique_code,
               reduction = reduction,
               tst = tst,
               ref = ref,
               grp = grp,
               ds = ds,
               covariates = covariates,
               test_use = test_use,
               assay = assay,
               table1 = table1,
               plot1 = plot1,
               plot2 = plot2)

}

#################################### Close log file
sink(type = "message")
sink()
close(logfile)