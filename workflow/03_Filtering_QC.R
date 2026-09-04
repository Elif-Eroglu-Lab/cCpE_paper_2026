############################## README
# Author: Inika Prasad
# Aim: To merge together spatial transcriptomics samples
# into the same Seurat object

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
# input_rds <- "../results/01_merge_samples/01_merged.rds"

# nFeature_threshold <- 200
# nCount_threshold <- 20000
# min_cells <- 3
# nsections <- 6

############################### Load Snakemake variables
input_rds <- snakemake@input[["input_rds"]]

nFeature_threshold <- snakemake@params[["nFeature_threshold"]]
nCount_threshold <- snakemake@params[["nCount_threshold"]]
min_cells <- snakemake@params[["min_cells"]]
nsections <- snakemake@params[["nsections"]]

plot0 <- snakemake@output[["plot0"]]
plot0_1 <- snakemake@output[["plot0_1"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
table <- snakemake@output[["table"]]
output_rds <- snakemake@output[["rds"]]


############################### Load data
obj <- readRDS(input_rds)

############################### Visualize cutoffs

# Visualize where the cutoffs will be
## Plot nFeature vs. nCount in cells, split by condition
nCount_vs_nFeature <- FeatureScatter(obj,
                                     shuffle = TRUE,
                                     feature1 = "nCount_Spatial",
                                     feature2 = "nFeature_Spatial",
                                     split.by = "condition") +
  geom_hline(yintercept = nFeature_threshold,
             linetype = "dashed",
             color = "darkgrey",
             linewidth = 1) +
  geom_vline(xintercept = nCount_threshold,
             linetype = "dashed",
             color = "darkgrey",
             linewidth = 1)

nCount_vs_nFeature
ggsave(filename = plot0, bg = "white")

## Plot nFeature vs. nCount in cells, split by batch
nCount_vs_nFeature <- FeatureScatter(obj,
                                     shuffle = TRUE,
                                     feature1 = "nCount_Spatial",
                                     feature2 = "nFeature_Spatial",
                                     split.by = "batch") +
  geom_hline(yintercept = nFeature_threshold,
             linetype = "dashed",
             color = "darkgrey",
             linewidth = 1) +
  geom_vline(xintercept = nCount_threshold,
             linetype = "dashed",
             color = "darkgrey",
             linewidth = 1)

nCount_vs_nFeature
ggsave(filename = plot0_1, bg = "white")

############################## Identify genes & cells to filter out

# Identify cells to filter out
selected_cells <- WhichCells(obj,
                             expression = nFeature_Spatial > nFeature_threshold & 
                                          nCount_Spatial < nCount_threshold)

# Identify genes to filter out:
countlayers <- obj@assays$Spatial@layers %>% names() #for every layer, i.e. sample

# Get genes that pass the cutoff for each sample
selft_per_sample <- list()
for (i in seq_along(countlayers)) {
  # Create a binary matrix of which cells express a gene at all
  binary_su <- obj[["Spatial"]][countlayers[i]] > 0

  # Identify features which pass the minimum cell number filter
  selected_fts <- rownames(binary_su)[Matrix::rowSums(binary_su) >= min_cells] %>% unlist()

  selft_per_sample[[countlayers[i]]] <- selected_fts

}

# Now selft_per_sample contains all genes which pass the minimum spot threshold
# in each sample
# Tally how many samples a gene is present in (above min sp)
freq_enough_gene <- table(unlist(selft_per_sample))
selected_features <- freq_enough_gene[freq_enough_gene >= nsections] %>% names()



#############################  Filter genes and cells
obj <- subset(obj, features = selected_features, cells = selected_cells)
# warnings here can be ignored: https://github.com/satijalab/seurat/issues/8848
joined_obj <- JoinLayers(obj)
############################# QC plot & save

# Compute new values of nFeature and nCount after filtering
# Because the pre-provided nFeature and nCount metrics don't auto-update at this stage
joined_obj$filt_nCount_Spatial <- colSums(joined_obj@assays$Spatial$counts)
joined_obj$filt_nFeature_Spatial <- colSums(joined_obj@assays$Spatial$counts != 0)

# Compute nFeature and nCount
VlnPlot(joined_obj, features = c("filt_nFeature_Spatial", "filt_nCount_Spatial"),
        ncol = 2,
        alpha = 0.1,
        group.by = "sample") + theme(rect = element_rect(fill = "transparent"))

ggsave(filename = plot1)

## View nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(joined_obj,
                                     feature1 = "filt_nCount_Spatial",
                                     feature2 = "filt_nFeature_Spatial",
                                     split.by = "condition")
nCount_vs_nFeature
ggsave(filename = plot2, bg = "white")

## Are there any genes which dominate the nCounts per cell? Overall
obj_countmat <- LayerData(joined_obj, assay = "Spatial", layer = "counts")
obj_countmat@x <- obj_countmat@x / rep.int(colSums(obj_countmat), diff(obj_countmat@p)) * 100
most_expressed <- order(Matrix::rowSums(obj_countmat), decreasing = TRUE)[20:1]
gene_names <- rownames(obj_countmat)[most_expressed]

# Make plot: genes that dominate nCounts per cell
png(plot3)
par(mar = c(5, 10, 4, 2)) # increase margin to prevent labels getting cut off
boxplot(as.matrix(t(obj_countmat[most_expressed, ])),
  cex = 0.1, las = 1, xlab = "Percent counts per spot",
  col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
  names = gene_names
)
dev.off()

## QC metrics
obj_summary <- joined_obj[[]] %>%
  group_by(sample) %>%
  summarise(
    spots = n(),
    total_features = dim(joined_obj)[1],
    mean_nCount = mean(filt_nCount_Spatial, na.rm = TRUE),
    mean_nFeatures = mean(filt_nFeature_Spatial, na.rm = TRUE)
  )

write_xlsx(x = obj_summary, path = table)

############################### Save filtered object

## Save filtered RDS object
saveRDS(object = obj, file = output_rds)

############################### Close log file
sink(type = "message")
sink()
close(logfile)