# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: For a chosen clustering resolution,
# generate list and plot of markers
# using FindAllMarkers (straighforward),
# or FindConservedMarkers (which finds markers conserved across datasets)

# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################# Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(tibble)
library(readxl)

############################ Set seed for reproducibility
set.seed(100)

############################ Load snakemake variables
rds_input <- snakemake@input[["rds"]]

table1 <- snakemake@output[["table1"]] # Excel of all markers
plot0 <- snakemake@output[["plot0"]] # UMAP showing clustering
plot1 <- snakemake@output[["plot1"]] # Dotplots

myreduction <- snakemake@params$myreduction # UMAP to be used for plots
myclusters <- snakemake@params$myclusters # Clustering that is being QC'd
grouping.var <- snakemake@params$grouping_var # metadata parameter for FindConservedMarkers

calc_conserved_markers <- snakemake@params$calc_conserved_markers

############################# Load data
obj <- readRDS(rds_input)

obj <- JoinLayers(obj)

source("workflow/functions.R")

################################# Check that inputs are valid
# Check that the chosen reduction exists
if (!myreduction %in% names(obj@reductions)) {
  stop(paste0(
    myreduction,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(obj@reductions), collapse = ", ")
  ))
}

# Check that the chosen clustering exists
if (!myclusters %in% colnames(obj@meta.data)) {
  stop(paste0(
    myclusters,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  ))
}

# Check that the chosen clustering exists
if (!grouping.var %in% colnames(obj@meta.data)) {
  stop(paste0(
    grouping.var,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  ))
}

############################# Print resolution
# Assign identities to cells based on the present cluster
Idents(obj) <- myclusters

# What levels does the clustering have?
clusters <- levels(obj)

# How many clusters do we have?
cli <- clusters %>% length()
cli

############################# Print UMAP
# So we know what clustering we're using
png(plot0)
DimPlot(obj, label = TRUE, reduction = myreduction)
dev.off()

############################# FindAllMarkers: WILCOX
# Find markers for every cluster using the Wilcox test
wilcox_sig <- FindAllMarkers(obj,
                                    test.use = "wilcox",
                                    only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25,
                                    max.cells.per.ident = 10000,
                                    min.cells.group = 10,
                                    random.seed = 100) %>%
  # filter only statistically significant markers
  filter(p_val_adj < 0.05) %>%
  ranker() %>%
  arrange(cluster, -score)

# Extract top 5 markers per cluster for dotplot
markers2plot <- wilcox_sig %>%
  group_by(cluster) %>%
  arrange(cluster, -score, , -(abs(avg_log2FC))) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

############################# FindConservedMarkers IF NEEDED

if (calc_conserved_markers == "TRUE") {

  ############################ CONSERVED MARKERS: ROC
  # Aim to find marker genes for queryclust regardless of grouping.var
  # Initialize list to store results
  conserved_markers_list_roc <- list()

  # Loop through clusters
  for (cluster in clusters) {
    message("Running FindConservedMarkers ROC on cluster: ", cluster)

    # Run FindConservedMarkers for each cluster
    conserved_markers <- FindConservedMarkers(
      obj,
      ident.1 = cluster,
      grouping.var = grouping.var,
      test.use = "roc",
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      max.cells.per.ident = 10000,
      min.cells.group = 10,
      random.seed = 100
    )

    # Store in list
    conserved_markers_list_roc[[cluster]] <- conserved_markers
  }

  # Combine all into one data frame with cluster label
  roc_all_conserved <- bind_rows(
    lapply(names(conserved_markers_list_roc), function(cluster) {
      df <- conserved_markers_list_roc[[cluster]]
      df <- rownames_to_column(df, var = "gene")
      df$cluster <- cluster
      df
    })
  )

  # Combine prediction power for both groups
  roc_all_conserved$total_power <- rowSums(roc_all_conserved[grep("power$", names(roc_all_conserved))])
  roc_all_conserved$pct_1all <- rowSums(roc_all_conserved[grep("pct.1$", names(roc_all_conserved))])
  roc_all_conserved$pct_2all <- rowSums(roc_all_conserved[grep("pct.2$", names(roc_all_conserved))])
  roc_all_conserved$abundance_diff <- roc_all_conserved$pct_1all - roc_all_conserved$pct_2all

  # Rank markers by importance
  roc_all_conserved <- roc_all_conserved %>%
    arrange(cluster, -total_power)

  # Subset markers for the dotplot
  markers3plot <- roc_all_conserved %>%
    group_by(cluster) %>%
    top_n(5, total_power) %>%
    arrange(cluster, -total_power) %>%
    pull(gene) %>%
    unique()

  ############################### CONSERVED MARKERS: WILCOX
  # Initialize list to store results
  conserved_markers_list_wilcox <- list()

  # Loop through clusters
  for (cluster in clusters) {
    message("Running FindConservedMarkers WILCOX on cluster: ", cluster)

    # Run FindConservedMarkers for each cluster
    conserved_markers <- FindConservedMarkers(
      obj,
      ident.1 = cluster,
      grouping.var = grouping.var,
      test.use = "wilcox",
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      max.cells.per.ident = 10000,
      min.cells.group = 10,
      random.seed = 100
    )

    # Store in list
    conserved_markers_list_wilcox[[cluster]] <- conserved_markers
  }

  # Combine all into one data frame with cluster label
  wilcox_all_conserved <- bind_rows(
    lapply(names(conserved_markers_list_wilcox), function(cluster) {
      df <- conserved_markers_list_wilcox[[cluster]]
      df <- rownames_to_column(df, var = "gene")
      df$cluster <- cluster
      df
    })
  )

  # Identify columns containing the adjusted p-values
  adj_pval.codes <- colnames(wilcox_all_conserved)[grepl(pattern = "*_p_val_adj$", x = colnames(wilcox_all_conserved))]

  # Compute the max adjusted p value from the given adj p-values
  wilcox_all_conserved$max_adj_pval <- apply(
    X = wilcox_all_conserved[, adj_pval.codes, drop = FALSE],
    MARGIN = 1, # in a given row
    FUN = max # maximum function
  )

  wilcox_all_conserved$pct_1all <- rowSums(wilcox_all_conserved[grep("pct.1$", names(wilcox_all_conserved))])
  wilcox_all_conserved$pct_2all <- rowSums(wilcox_all_conserved[grep("pct.2$", names(wilcox_all_conserved))])
  wilcox_all_conserved$abundance_diff <- wilcox_all_conserved$pct_1all - wilcox_all_conserved$pct_2all

  # Subset the significant genes
  wilcox_sig_conserved <- wilcox_all_conserved %>%
    # filter only statistically significant markers
    filter(max_adj_pval < 0.05) %>%
    arrange(cluster)

  markers4plot <- wilcox_sig_conserved %>%
    group_by(cluster) %>%
    top_n(5, abundance_diff) %>%
    arrange(cluster, -abundance_diff) %>%
    pull(gene) %>%
    unique()

  # Save the markers
  writexl::write_xlsx(list(wilcox_sig = wilcox_sig,
                           roc_all_conserved = roc_all_conserved,
                           wilcox_all_conserved = wilcox_all_conserved,
                           wilcox_sig_conserved = wilcox_sig_conserved),
                      table1)

  setsofmarkers <- list(wilcox_sig = markers2plot,
                        roc_all_conserved = markers3plot,
                        wilcox_sig_conserved = markers4plot)

} else {

  # Save marker lists
  writexl::write_xlsx(list(wilcox_sig = wilcox_sig),
                      table1)

  setsofmarkers <- list(wilcox_sig = markers2plot)

}


################################ Plotting significant markers

## Print dotplot
pdf(plot1,
    height = 3 + (as.integer(cli) * 0.6),
    width  = 4 + (as.integer(cli) * 0.4))

for (markerset in setsofmarkers) {
  p <- DotPlot(obj,
               features = markerset,
               cols = "RdBu") +
    coord_flip() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      axis.title = element_blank()
    )

  print(p)

}

dev.off()

################################ Close log file
# Close log file
sink(type = "message")
sink()
close(logfile)