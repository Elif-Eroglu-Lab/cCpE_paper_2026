# README
# Author: Inika Prasad, inika.prasad@gmail.com
# Aim: For a chosen clustering resolution,
# See how the sample split looks across clusters,
# and how many doublets we get for clusters,
# with different visualizations
# And check how the known markers look across clusters
# using dotplots

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
library(writexl)
library(viridis)

set.seed(100)

################################# Tester variables
myreduction <- "umap.pca"
myclusters <- "pca_resolution_1.2"
knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

################################# Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot1_1 <- snakemake@output[["plot1_1"]]
plot1_2 <- snakemake@output[["plot1_2"]]
plot3 <- snakemake@output[["plot3"]]
plot4 <- snakemake@output[["plot4"]]
table1 <- snakemake@output[["table1"]]
plot5 <- snakemake@output[["plot5"]]
plot6 <- snakemake@output[["plot6"]]
plot7 <- snakemake@output[["plot7"]]
plot8 <- snakemake@output[["plot8"]]
plotx <- snakemake@output[["plotx"]]

myreduction <- snakemake@params$myreduction # UMAP to be used for plots
myclusters <- snakemake@params$myclusters # Clustering that is being QC'd
knownmarker_file <- snakemake@params$knownmarker_file
knownmarker_id_col <- snakemake@params$knownmarker_id_col

################################# Load dataset
obj <- readRDS(input_rds)

################################# Check that inputs are valid
# Check that the chosen reduction and clustering exist,
# Abort with a message if they don't

if (!myreduction %in% names(obj@reductions)) {
  stop(paste0(
    myreduction,
    " is not present in the Seurat object.\n",
    "Available reductions are: ",
    paste(names(obj@reductions), collapse = ", ")
  ))
}

if (!myclusters %in% colnames(obj@meta.data)) {
  stop(paste0(
    myclusters,
    " is not present in the Seurat object.\n",
    "Available metadata/clusterings are: ",
    paste(colnames(obj@meta.data), collapse = ", ")
  ))
}

if (!file.exists(knownmarker_file)) {
  stop(paste0(
    knownmarker_file, " does not exist. Please check path."
  ))
}

################################# Dotplot known marker genes
# Specify genes of interest (known marker genes)
goi <- read.csv(knownmarker_file, sep = ";")
print(goi)

# Clean up list
goi <- goi[goi[[knownmarker_id_col]] != "",]
goi <- goi[!duplicated(goi[[knownmarker_id_col]]), ]

# Create list of named lists
celltype_lists <- goi %>%
  group_by(Cell.Type) %>%
  summarise(gene = list(!!sym(knownmarker_id_col)), .groups = "drop") %>%
  { setNames(.$gene, .$Cell.Type) }

# Reorder list
categories_order <- unique(goi$Cell.Type)

celltype_lists <- setNames(
  lapply(categories_order, function(ct) goi[[knownmarker_id_col]][goi$Cell.Type == ct]),
  categories_order
)

pdf(plot0, width = 30, height = 8.3)

# Scaled plot
DotPlot(obj,
        group.by = myclusters,
        features = celltype_lists,
        cols = "RdBu",
        scale = TRUE) +
  RotatedAxis()

# Unscaled plot
DotPlot(obj,
        group.by = myclusters,
        features = celltype_lists,
        cols = "RdBu",
        scale = FALSE) +
  RotatedAxis()

dev.off()

################################## QC plots

# Plot sample of origin
p1 <- DimPlot(obj,
              reduction = myreduction,
              group.by = "orig.ident",
              label.size = 2) +
  ggtitle("Samples") +
  theme(legend.position = "bottom")
p1
ggsave(plot1)

# Plot clusters
p2 <- DimPlot(obj,
              reduction = myreduction,
              group.by = myclusters,
              label = TRUE) +
  ggtitle(myclusters)
p2
ggsave(plot2, width = 8)

# Split into samples, monochrome
p1_1 <- DimPlot(obj,
                reduction = myreduction,
                split.by = "orig.ident",
                group.by = "orig.ident") +
  ggtitle("Samples") +
  theme(legend.position = "bottom")
p1_1
ggsave(plot1_1, width = 13)

# Split into samples, color by clusters
p1_2 <- DimPlot(obj,
                reduction = myreduction,
                split.by = "orig.ident",
                group.by = myclusters,
                label = TRUE) +
  ggtitle(myreduction)
p1_2
ggsave(plot1_2, width = 18)

# Plot doublet confidence
p3 <- DimPlot(obj,
              reduction = myreduction,
              group.by = "Doublet_Confidence") +
  theme(legend.position = "bottom")
p3
ggsave(plot3)

# Plot nCount and nFeature
p4 <- FeaturePlot(obj,
                  reduction = myreduction,
                  features = c("nCount_RNA", "nFeature_RNA"),
                  label.size = 2)
p4
ggsave(plot4, width = 16)

# Plot Cell cycle scores
p7 <- FeaturePlot(obj,
                  reduction = myreduction,
                  features = c("S.Score", "G2M.Score"),
                  label.size = 2)
p7
ggsave(plot7, width = 16)

# Plot Cell cycle stage categorization
p8 <- DimPlot(obj,
              reduction = myreduction,
              group.by = "Phase") +
  theme(legend.position = "bottom")
p8
ggsave(plot8)


############################################ QC: no. of cells across clusters & samples
# No. of cells per sample per cluster
df <- table(obj@meta.data[["orig.ident"]], obj@meta.data[[myclusters]])

# Each sample is a row. Remove samples that have no cells at all
# Since I defined the levels for all samples earlier, and obj
# may not contain cells from all the samples
df <- df[rowSums(df) != 0, ] # first row

df2 <- df %>%
  as.data.frame.matrix() %>%
  rownames_to_column(., var = "sample")

############################################ QC Plot: sample proportions across clusters
# all clusters in a sample add up to 100% (account for different no. of cells per sample)
propmat <- prop.table(df, margin = 1)

df3 <- propmat %>%
  as.data.frame.matrix() %>%
  rownames_to_column(., var = "sample")

# proportions per cluster add up to 100% (account for cluster size)
propmat <- prop.table(propmat, margin = 2) * 100

df4 <- propmat %>%
  as.data.frame.matrix() %>%
  rownames_to_column(., var = "sample")

# Save no. of cells per sample per cluster
write_xlsx(list(ncells = df2,
                ncells_normrow = df3,
                ncells_normrowcol = df4),
           path = table1)

# Prep data for heatmap
df2plot <- propmat %>% as.data.frame()
colnames(df2plot) <- c("Sample", "Cluster", "Proportion")

######################################
# Make a barplot
df3_long <- df3 %>%
  pivot_longer(!sample, names_to = "cluster", values_to = "percentage")

# raw values (ncells)
df2_long <- df2 %>%
  pivot_longer(!sample, names_to = "cluster", values_to = "cellnumber")

df_plot <- left_join(df3_long, df2_long, by = c("cluster", "sample"))


# df_plot$sample <- factor(df_plot$sample, levels = c("10dpci_veh", "10dpci_cCpE"))

# Plot total proportions
ggplot(df_plot, aes(fill = sample, y = percentage, x = cluster)) +
  geom_bar(position = position_dodge(width = 0.6), width=0.6, stat = "identity") +
  ylab("Percentage of total population") +
  xlab("Cluster") +
  coord_flip() +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  geom_text(aes(label = cellnumber),
            position = position_dodge(width = 0.9),
            hjust = -0.2,       # adjust depending on preference
            size = 3)
ggsave(plotx)

#######################################
# Get numbers and proportions of cells per cluster
print(df)
print(df2plot)

# hierachical clustering for heatmap
dist_cols <- dist(t(propmat)) # compute distance
hc_cols <- hclust(dist_cols) # cluster columns based on distance
category_order <- hc_cols$labels[hc_cols$order] #reorder labels based on clustering
df2plot$Cluster <- factor(df2plot$Cluster, levels = category_order) #reorder Cluster levels

# Heatmap of proportions
rd_bu_colors <- rev(brewer.pal(5, "RdBu")) # choosing 3 makes the heatmap paler
ggplot(df2plot, aes(x = Cluster, y = Sample, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = rd_bu_colors,
    values = scales::rescale(c(0, 25, 100)),  # midpoint at 0.25
    limits = c(0, 100),
    name = "Proportion (%)"
  ) +
  theme_minimal() +
  xlab(myclusters) +
  ggtitle("Sample proportions per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

nclust <- df2plot$Cluster %>% unique() %>% length()
ggsave(plot5, width = 6 + (nclust * 0.14), height = 2, bg = "white")

############################################ QC Plot: doublets per cluster
# proportion predicted doublets per cluster
df <- table(obj@meta.data[["Doublet_Confidence"]], obj@meta.data[[myclusters]])
propmat <- prop.table(df, margin = 2) * 100
df2plot <- propmat %>% as.data.frame()
colnames(df2plot) <- c("Doublet_Confidence", "Cluster", "Proportion")

# hierachical clustering for heatmap
dist_cols <- dist(t(propmat)) # compute distance
hc_cols <- hclust(dist_cols) # cluster columns based on distance
category_order <- hc_cols$labels[hc_cols$order] #reorder labels based on clustering
df2plot$Cluster <- factor(df2plot$Cluster, levels = category_order) #reorder Cluster levels

# Heatmap of proportions
rd_bu_colors <- rev(brewer.pal(5, "RdBu"))
ggplot(df2plot, aes(x = Cluster, y = Doublet_Confidence, fill = Proportion)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colors = rd_bu_colors,
    values = scales::rescale(c(0, 50, 100)),
    limits = c(0, 100),
    name = "Proportion (%)"
  ) +
  theme_minimal() +
  xlab(myclusters) +
  ggtitle("Doublet proportions per Cluster") +
  theme(plot.title = element_text(hjust = 0.5))

nclust <- df2plot$Cluster %>% unique() %>% length()
ggsave(plot6, width = 6 + (nclust * 0.14), height = 2, bg = "white")

#################################### Close log file
sink(type = "message")
sink()
close(logfile)