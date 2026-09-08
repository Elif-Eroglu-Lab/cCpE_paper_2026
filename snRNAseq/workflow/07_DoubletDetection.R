# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)
library(tibble)
library(tidyr)

if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
}

library(DoubletFinder)

set.seed(100)

# Load data
sfi <- readRDS(snakemake@input[["rds"]])

# Pull metadata
sample <- snakemake@params[["sample"]]
ndims <- as.numeric(snakemake@params[["ndims"]])
chosen_resolution <- snakemake@params[["chosen_resolution"]]

print(paste("ndims for", sample, "is", ndims, sep = " "))

# Calculate how BC metrics change with different parameters
sweep.res <- paramSweep(sfi, PCs = 1:ndims)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Save plot: BC metrics against the pK parameter IDs
old_par <- par(no.readonly = TRUE)
par(mgp = c(3.5, 1, 0)) # temporarily increase distance between y axis & ticks
png(filename = snakemake@output[["plot1"]])
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las = 2,
        xlab = "pK", ylab = "BCmetric", main = sample,
        ylim = c(0, max(bcmvn$BCmetric)), axis.lty = 1)
dev.off()
par(old_par)

# Choose the pK that is associated with the maximum BC metric
pK_selected <- as.vector(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)] %>% as.numeric
print(paste("pK_selected for", sample, "is", pK_selected, sep = " "))

# define the expected number of doublets.
nExp <- round(ncol(sfi) * snakemake@params[["doublet_rate"]]) # expected % doublets
print(paste("nExp for", sample, "is", nExp, sep = " "))

# Find the doublets (without correcting for homotypic doublets)
sfi_mr <- doubletFinder(sfi, pN = 0.25, pK = pK_selected, nExp = nExp, PCs = 1:ndims)
print("Completed doubletFinder without homotypic doublet proportion adjustment")

# name of the DF prediction & pANN can change, so extract the correct column name
DF.name <- colnames(sfi_mr@meta.data)[grepl("DF.classification", colnames(sfi_mr@meta.data))]
print(paste("DF.name for", sample, "is", DF.name, sep = " "))

pANN.name <- colnames(sfi_mr@meta.data)[grepl("pANN", colnames(sfi_mr@meta.data))]
print(paste("pANN.name for", sample, "is", pANN.name, sep = " "))

## Estimate proportion of homotypic doublets
annotations <- sfi@meta.data[[chosen_resolution]]
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(snakemake@params[["doublet_rate"]] * nrow(sfi@meta.data))
nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
print(paste("nExp_poi.adj for", sample, "is", nExp_poi.adj, sep = " "))

# Find the doublets after correcting for homotypic doublets
sfi_mr <- doubletFinder(sfi_mr, pN = 0.25, pK = pK_selected, nExp = nExp_poi.adj, PCs = 1:ndims, reuse.pANN = pANN.name)
print("Completed doubletFinder with homotypic doublet proportion adjustment")

colnames(sfi_mr@meta.data)

# Extract the column name for homotypic doublet adjusted classification
DF.names <- colnames(sfi_mr@meta.data)[grepl("DF.classification", colnames(sfi_mr@meta.data))]
DF.names
adjusted_DF.name <- DF.names[DF.names != DF.name]
print(paste("adjusted_DF.name for", sample, "is", adjusted_DF.name, sep = " "))

# Create metadata for confidence of doublet classification
# Column names made evaluable
colA_sym <- sym(DF.name)
colB_sym <- sym(adjusted_DF.name)

sfi_mr@meta.data <- sfi_mr@meta.data %>%
  mutate(Doublet_Confidence = case_when(
    !!colA_sym == "Doublet" & !!colB_sym == "Doublet" ~ "High_Confidence_Doublet",
    xor(!!colA_sym == "Doublet", !!colB_sym == "Doublet") ~ "Low_Confidence_Doublet",
    !!colA_sym == "Singlet" & !!colB_sym == "Singlet" ~ "Singlet",
    TRUE ~ NA_character_
  ))

# Save plot: UMAP with the doublets highlighted
png(filename = snakemake@output[["plot2"]], width = 600)
DimPlot(sfi_mr, group.by = "Doublet_Confidence") +
ggtitle(paste0(sample, ": Doublet Confidence")) +
theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Save plot: Violin with cells & their no. of Features, split by singlet/doublet
png(filename = snakemake@output[["plot3"]], width = 600)
VlnPlot(sfi_mr, features = c("nFeature_RNA", "nCount_RNA"),
        group.by = "Doublet_Confidence", pt.size = 0.1, alpha = 0.2)
dev.off()

# Save plot: pANN score across clusters
png(filename = snakemake@output[["plot4"]])
VlnPlot(sfi_mr, features = pANN.name[1], group.by = chosen_resolution, alpha = 0.2)
dev.off()

# Confirm cluster identitities, ease visualization
png(filename = snakemake@output[["plot5"]])
DimPlot(sfi_mr, group.by = chosen_resolution, label = TRUE) +
ggtitle(sample) +
theme(plot.title = element_text(hjust = 0.5))
dev.off()

## View nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(sfi_mr,
                                     feature1 = "nCount_RNA",
                                     feature2 = "nFeature_RNA",
                                     group.by = "Doublet_Confidence",
                                     shuffle = TRUE) +
  ggtitle(snakemake@wildcards[["sample"]])

ggsave(filename = snakemake@output[["plot6"]], bg = "white")

# Extract numbers of doublets & singlets
x <- table(sfi_mr[["Doublet_Confidence"]])
x

# Save the details of the doublet finding & results
results_df <- data.frame(sample = sample,
                         clustering = chosen_resolution,
                         pK = pK_selected,
                         nExp = nExp,
                         nExp_poi.adj = nExp_poi.adj,
                         DF.name = DF.name,
                         homotypic_adjusted_DF.name = adjusted_DF.name,
                         pANN.name = pANN.name[1],
                         predicted_doublets_high_confidence = x["High_Confidence_Doublet"], 
                         predicted_doublets_low_confidence = x["Low_Confidence_Doublet"],
                         predicted_singlets = x["Singlet"])

write_xlsx(x = results_df, path = snakemake@output[["table"]])

# Calculate percentage contribution of doublets to each cluster
df <- table(sfi_mr[["Doublet_Confidence", drop = TRUE]], sfi_mr[[chosen_resolution, drop = TRUE]]) %>%
  as.data.frame.matrix() %>%
  mutate(across(everything(), ~ . / sum(.) * 100)) %>%
  rownames_to_column(var = "Proportions")

# Put into long format so that ggplot can plot
long_spots <- pivot_longer(df,
                           cols = -Proportions,
                           names_to = "Cluster",
                           values_to = "pct_doublets")

long_spots$Cluster <- factor(long_spots$Cluster, levels = unique(long_spots$Cluster) )
levels(long_spots$Cluster)

# Make plot for proportions
p <- ggplot(long_spots, aes(x = Cluster, y = pct_doublets, fill = Proportions))+
  geom_bar(stat = "identity", colour = "black") +
  ylab("Percentage") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal() +
  theme(
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  ggtitle(label = sample,
          subtitle = "Proportions of doublets per cluster") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p

ggsave(filename = snakemake@output[["plot7"]], bg = "white")

# Save the Seurat object with doublet/singlet annotations
saveRDS(sfi_mr, snakemake@output[["rds"]])

# Close log file
sink(type = "message")
sink()
close(logfile)