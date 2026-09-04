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

source("workflow/functions.R")

################################# Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

unique_code <- snakemake@params$unique_code
input_rds <- snakemake@params$input_rds
clustering <- snakemake@params$clustering
reduction <- snakemake@params$reduction
tst <- snakemake@params$tst
ref <- snakemake@params$ref
grp <- snakemake@params$grp
ds <- snakemake@params$ds
covariates <- snakemake@params$covariates
test_use <- snakemake@params$test_use
assay <- snakemake@params$assay

# Outputs
table1 <- snakemake@output[["table1"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]

################################ Define hardcoded variables
input_rds <- snakemake@input[["rds"]]

################################# Load dataset
obj <- readRDS(input_rds)

# Need to ensure this is already done.
obj2 <- JoinLayers(obj)

# Source the file where the metadata columns are created
source("workflow/08b_make_metadata4DEGs.R")

################################# Check Snakemake variables

# Check if grouping exists in Seurat object
if (!grp %in% colnames(obj2@meta.data)) {
  stop(paste0(
    grp,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj2@meta.data), collapse = ", ")
  ))
}

available_levels <- obj2@meta.data[[grp]] %>% unique()
levels2check <- c(ref, tst)

# Check if ref & tst are present in that grouping
if (!all(levels2check %in% available_levels)) {
  missing_levels <- levels2check[!levels2check %in% available_levels]
  stop(paste0(
    "The following ref/tst is not present in the metadata column called ",
    grp,
    ": ",
    paste(missing_levels, collapse = ", "), "\n",
    "Available ref/tst are: ",
    paste(available_levels, collapse = ", ")
  ))
}

# Parse covariates provided
if (!is.null(covariates)) {
  latent.vars <- covariates %>% strsplit(split = ",") %>% unlist() %>% trimws()

    # Check if covariates are present in Seurat object
  if (!all(latent.vars %in% colnames(obj2@meta.data))) {
    missing_vars <- latent.vars[!latent.vars %in% colnames(obj2@meta.data)]
    stop(paste0(
      "The following covariates are not present in the metadata column:",
      paste(missing_vars, collapse = ", "), "\n",
      "Available metadata/clusterings are: ",
      paste(colnames(obj2@meta.data), collapse = ", ")
    ))
  }
} else {
  latent.vars <- NULL
}

################################ DGE from scratch

description <- paste0("DGE-", grp, "-", ref, "_vs_", tst, "_", covariates)
print(description)

# Copy of Seurat object to subset & use for DGE
Idents(obj2) <- grp
obj_dge <- obj2

# Assign identities to Seurat object based on grouping
Idents(obj_dge) <- grp

# check if ref and test levels are present in grp
Idents(obj_dge) %>% unique()

# Subset object to include test and reference cells
obj_dge <- subset(obj_dge, subset = !!sym(grp) %in% c(ref, tst))

# Save this non-downsampled original object for subsequent heatmap generation
obj_nods <- obj_dge

# Find number to downsample to, if downsampling 200 is not possible
ncells_per_grp <- Idents(obj_dge) %>% table %>% as.data.frame()
colnames(ncells_per_grp) <- c("group", "cells")
ncells_per_grp


if (any(ncells_per_grp$cells < ds)) {
  new_ds <- min(ncells_per_grp$cells)
  message(paste0("One or more groups have cells less than ",
                 ds,
                 " cells. Downsampling to the minimum no. of cells instead, 
                 which is ",
                 new_ds))
  print(ncells_per_grp)

  ds <- new_ds

  # Do the downsampling
  message(paste0("Downsampling data to ", ds, " cells per group."))
  obj_dge <- subset(obj_dge, downsample = ds)

} else if(any(ncells_per_grp$cells >= ds)) {
  # Do the downsampling
  message(paste0("Downsampling data to ", ds, " cells per group."))
  obj_dge <- subset(obj_dge, downsample = ds)

} else if (is.null(ds)) {
  ds <- "no downsampling"
} else {
  stop("Error in downsampling")
}


# Get count data for both groups (after possible downsampling)
ncells_used <- Idents(obj_dge) %>% table %>% as.data.frame()
colnames(ncells_used) <- c("group", "cells")
ncells_used

# Must have at least 3 cells per group
if (sum(ncells_used$cells < 3) != 0) {
  stop("There are less than 3 cells in at least one group")
}

# Do the DGE testing
m <- FindMarkers(obj_dge,
                 ident.1 = tst,
                 ident.2 = ref,
                 group.by = grp,
                 test.use = test_use,
                 only.pos = FALSE,
                 min.pct = 0.25,
                 logfc.threshold = 0,
                 min.cells.group = 3,
                 random.seed = 100,
                 assay = assay,
                 slot = "data",
                 latent.vars = latent.vars) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, everything()) %>%
  ranker()

run_on <- Sys.time()

# Subset DEGs to include only significant adj.p.val < 0.05
sig_degs <- m %>% filter(p_val_adj < 0.05)

# Construct summary
dlist <- data.frame(
  unique_code = unique_code,
  grp = grp, ref = ref, test = tst,
  cells_ref = ncells_used[ncells_used$group == ref, ]$cells,
  cells_test = ncells_used[ncells_used$group == tst, ]$cells,
  degs_up = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC > 0)),
  degs_down = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC < 0)),
  degs_tot = nrow(dplyr::filter(m, p_val_adj < 0.05)),
  skewness = moments::skewness(m$p_val),
  kurtosis = moments::kurtosis(m$p_val),
  description = description,
  statistical_test = test_use,
  downsampling = ds,
  covariates = paste(latent.vars, collapse = ", "),
  run_on = run_on,
  comment = NA,
  clustering = clustering
)

# Save info
writexl::write_xlsx(list(full_degs = m,
                         sig_degs = sig_degs,
                         run_details = dlist),
                    table1)

# Top 10 genes for violin plots
m1 <- m %>%
  ranker() %>%
  slice_head(n = 10)

# Cells used for the DGE in each group
ref_cells <- colnames(obj_dge)[obj_dge[[grp]] == ref]
tst_cells <- colnames(obj_dge)[obj_dge[[grp]] == tst]

# Save result plots in a pdf
pdf(plot1, height = 4, width = 4)
plot_volcano(m, ref_grp = ref, test_grp = tst)

p1 <- DimPlot(obj2, cells.highlight = ref_cells, reduction = reduction) +
  ggtitle(ref) +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(obj2, cells.highlight = tst_cells, reduction = reduction) +
  ggtitle(tst) +
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5))

p1
p2

plot_pvaldist(m) + ggtitle("Raw p-value distribution")
plot_qq(m) + ggtitle("QQ plot")
plot_violin_seurat(obj_dge, grp, c(ref, tst),
                   g = m1$gene, assay = assay, ncol = 2, base_size = 6)
dev.off()

# Make a dotplot showing where the genes are expressed
pdf(plot2, height = 2.5 + (0.2 * nrow(sig_degs)), width = 4 + 0.4 * length(unique(Idents(obj2))))
DotPlot(obj2, features = sig_degs$gene, cols = "RdBu", scale = TRUE) + RotatedAxis() + coord_flip()
DotPlot(obj2, features = sig_degs$gene, cols = "RdBu", scale = FALSE) + RotatedAxis() + coord_flip()
dev.off()

#################################### Close log file
sink(type = "message")
sink()
close(logfile)