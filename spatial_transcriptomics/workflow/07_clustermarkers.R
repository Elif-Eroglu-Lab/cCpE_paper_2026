# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aims: To find markers defining each cluster &
# plot the top markers.

# User must define the name of the clustering column
# in the metadata
# e.f. cca_resolution_1

############################## Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

############################## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(pals)
library(patchwork)

# Set seed for reproducibility
set.seed(100)

# Load helper functions
source("workflow/functions.R")

############################## Load Snakemake variables
rds_input <- snakemake@input[["rds"]]
res <- snakemake@params[["resolution"]]

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
table_output <- snakemake@output[["table"]]
# sample <- snakemake@params[["sample"]]

############################### Load data
sfi <- readRDS(rds_input)

############################### Check basics of data & clustering
print(res)

# Assign identities to cells based on the present cluster
Idents(sfi) <- sfi@meta.data[[res]]

# How many clusters do we have?
cli <- levels(sfi@active.ident) %>% length()

# High-contrast color palette
# selected_cols <- pals::glasbey(30)
# names(selected_cols) <- Idents(sfi) %>% levels()

# High-contrast color palette
selected_cols <- pals::glasbey(30)
names(selected_cols) <- sfi@meta.data[[res]] %>% unique()

pdf(plot0)
# Print clusters on sections
p_spatial <- plotall(sfi,
                     selected_clustering = res,
                     selected_cols = selected_cols,
                     crop = TRUE,
                     pt.size = 3)

combined_spatial <- wrap_plots(p_spatial, ncol = 3, guides = "collect")
combined_spatial
dev.off()


## Print UMAP
png(plot1)
DimPlot(sfi, label = TRUE, cols = selected_cols)
dev.off()

## Print UMAP, split by condition
png(plot3, width = 960)
DimPlot(sfi, label = TRUE, split.by = "condition", cols = selected_cols)
dev.off()


############################### Find Markers

# Find markers for every cluster using the Wilcox test
markers <- FindAllMarkers(sfi,
                          test.use = "wilcox", only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, max.cells.per.ident = 200,
                          min.cells.group = 10, random.seed = 100) %>%
  # filter only statistically significant markers
  filter(p_val_adj < 0.05) %>%
  ranker() %>%
  arrange(cluster, -score)

# Add a score that accounts for relative abundance
# should be a large positive value for the most specific markers
markers$abundancediff <- markers$pct.1 - markers$pct.2

# Create a new score metric
markers$score2 <- markers$score * markers$abundancediff

# Group cells by cluster and extract top 5 markers
markers2plot <- markers %>%
  group_by(cluster) %>%
  top_n(5, score) %>%
  arrange(cluster, -score)

############################### Plot top markers
p <- DotPlot(sfi, features = unique(markers2plot$gene),
             cols = "RdBu") +
  coord_flip() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.title = element_blank()
  )

ggsave(plot2, p, height = 3 + (as.integer(cli) * 0.6), width = 4 + (as.integer(cli) * 0.4), bg = "white")

############################### Save marker list
write_xlsx(x = markers, path = table_output)

############################### Close log file
sink(type = "message")
sink()
close(logfile)