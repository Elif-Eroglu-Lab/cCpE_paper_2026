#################################### README
# Author: Inika Prasad
# Created on: 2025-11-24
# Aim: # Evaluate expression of single-nuclei top markers in query transcriptomics

# to see overlap between the top markers of different clusters
# and to see if the unmapped clusters from workflow/09a_deconvolution.R
#  have poorly expressed top markers

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

#################################### Define variables
# 10dpci epicardial peelout
reference_markers <- "../../snRNAseq_2025/code/snakemake/results/16b_clustermarkers_4peeledintdata/10dpci/10dpci_DRscrub_Epi_markers.xlsx"
sheet <- "wilcox_sig_conserved"

# 10dpci markers
reference_markers <- "../../snRNAseq_2025/code/snakemake/results/12b_clustermarkers_4intdata/10dpci_DRscrub_markers.xlsx"
sheet <- "wilcox_sig_conserved"

# query dataset
query_path <- "../results/06_dimred_clustering/dimred.rds"
query_path <- "../results/06_dimred_clustering/allsections/allsections_dimred.rds"

ntopgenes <- 10

#################################### Load data
# Load marker lists from reference
markers <- read_xlsx(reference_markers,
                     sheet = sheet)

# Load query object
query <- readRDS(query_path)

# Sanity check query object
DimPlot(query)

#################################### Identify markers present in query
markers$inquery <- markers$gene %in% rownames(query@assays$Spatial)

# Check if the no. of markers per cluster OR no. of markers available in query is the cause
markers %>%
  group_by(cluster) %>%
  summarise(n = n(),
            inquery_sum = sum(inquery, na.rm = TRUE)) %>%
  View()


####################################### Extract top markers for module scores
cluster_gene_vectors <- markers %>%
  filter(inquery == TRUE) %>%
  group_by(cluster) %>%
  slice_max(order_by = abundance_diff, n = ntopgenes, with_ties = FALSE) %>%
  arrange(cluster, desc(abundance_diff)) %>%
  summarise(genes = list(gene), .groups = "drop") %>%
  deframe()

cluster_gene_vectors

##################################### Add Module Scores to query
query <- AddModuleScore(query,
                        features = cluster_gene_vectors,
                        name = "module")

# identify the metadata columns, which all starr with "name"
module_names <- colnames(query@meta.data)[grep("module", colnames(query@meta.data))]
module_names
##################################### Plot module scores globally
# Per condition
VlnPlot(query,
        features = module_names,
        group.by = "condition",
        alpha = 0.01)

# Per sample
VlnPlot(query,
        features = module_names,
        group.by = "sample",
        assay = "Spatial",
        alpha = 0.01)

#################################### Plot module scores in UMAP: all

FeaturePlot(query, features = module_names, cols = c("grey", "darkred"))

#################################### Plot module scores in UMAP: top

# Plot the top x percentile of cells per module score

for (module in module_names) {
  print(module)
  threshold <- quantile(query@meta.data[[module]], 0.90, na.rm = TRUE)
  print(threshold)
  topQuartile <- paste0(module, "_topQuartile")
  print(topQuartile)
  query@meta.data[[topQuartile]] <- query@meta.data[[module]] >= threshold
}

# Find metadata columns with the top quartile
moduletopQ_names <- colnames(query@meta.data)[grep("^module.*topQuartile$", colnames(query@meta.data))]

DimPlot(query, group.by = moduletopQ_names, split.by = "condition")

################################## Compare conditions
# Compare how many top percentile spots are in the different conditions

# Get info about top scorers and conditions
mods_table <- query@meta.data[c("condition", moduletopQ_names)]

# Pivot to a longer form of the table
df_long <- mods_table %>%
  pivot_longer(
    cols = moduletopQ_names,
    names_to = "moduletopQ",
    values_to = "value"
  )

# Make stacked barchart
ggplot(df_long, aes(x = condition, fill = value)) +
  geom_bar(position = "fill") +
  facet_wrap(~ moduletopQ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Condition",
    y = "Percent",
    fill = "Spots scoring 90% in module",
    title = "Module scoring by Condition"
  ) +
  theme_minimal(base_size = 14)

#################################### Plot module scores spatially (if query is spatial)
x <- plotall(query,
        crop = TRUE,
        features2plot = module_names,
        image.alpha = 0,
        pt.size.factor = 2)
x[[14]]

#####################################