###################################### README
# Aim: to filter out cells with too few features,
# and genes that are not abundantly enough expressed

####################################### Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

####################################### Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

###################################### Load Snakemake variables
input_rds <- snakemake@input[["rds"]]

nFeature_threshold <- snakemake@params[["nFeature_threshold"]]
nCount_threshold <- snakemake@params[["nCount_threshold"]]
min_cells <- snakemake@params[["min_cells"]]
sample <- snakemake@wildcards[["sample"]]

plot0 <- snakemake@output[["plot0"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
table <- snakemake@output[["table"]]
output_rds <- snakemake@output[["rds"]]

####################################### Load data
su <- readRDS(input_rds)

####################################### Define cutoffs
## Plot nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(su,
                                     feature1 = "nCount_RNA",
                                     feature2 = "nFeature_RNA") +
  ggtitle(sample) +
  geom_hline(yintercept = nFeature_threshold, linetype="dashed", color = "grey") +
  geom_vline(xintercept = nCount_threshold)

ggsave(filename = plot0, bg = "transparent")


####################################### Identify cells & features to filter
binary_su <- su[["RNA"]]$counts > 0

selected_cells <- WhichCells(su, expression = nFeature_RNA > nFeature_threshold & 
                                              nCount_RNA < nCount_threshold)
selected_features <- rownames(su)[Matrix::rowSums(binary_su) >= min_cells]

####################################### Filter cells and genes

sf <- subset(su, features = selected_features, cells = selected_cells)

####################################### Update the nCount & nFeatures
# Since they don't update automatically after subsetting
sf$filt_nCount_RNA <- colSums(sf@assays$RNA$counts)
sf$filt_nFeature_RNA <- colSums(sf@assays$RNA$counts != 0)

####################################### QC plot & save
VlnPlot(sf, features = c("filt_nFeature_RNA", "filt_nCount_RNA"), ncol = 2, alpha = 0.1)
ggsave(filename = plot1, bg = "transparent")

## View nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(sf,
                                     feature1 = "filt_nCount_RNA",
                                     feature2 = "filt_nFeature_RNA") +
                                     ggtitle(sample)

ggsave(filename = plot2, bg = "transparent")

## Are there any genes which dominate the nCounts per cell?
sf_countmat <- LayerData(su, assay = "RNA", layer = "counts")
sf_countmat@x <- sf_countmat@x / rep.int(colSums(sf_countmat), diff(sf_countmat@p)) * 100
most_expressed <- order(Matrix::rowSums(sf_countmat), decreasing = TRUE)[20:1]
gene_names <- rownames(sf_countmat)[most_expressed]

# Make plot: genes that dominate nCounts per cell
png(plot3)
par(mar = c(5, 10, 4, 2)) # increase margin to prevent labels getting cut off
boxplot(as.matrix(t(sf_countmat[most_expressed, ])),
  cex = 0.1, las = 1, xlab = "Percent counts per cell",
  col = (scales::hue_pal())(20)[20:1], horizontal = TRUE,
  names = gene_names
)
dev.off()

## QC metrics
sf_summary <- sf[[]] %>%
  summarise(
    cells = n(),
    total_features = dim(sf)[1],
    mean_nCount = mean(filt_nCount_RNA, na.rm = TRUE),
    mean_nFeatures = mean(filt_nFeature_RNA, na.rm = TRUE)
  )

###################################### Save QC metrics
write_xlsx(x = sf_summary, path = table)

###################################### Save RDS object
saveRDS(object = sf, file = output_rds)

###################################### Close log file
sink(type = "message")
sink()
close(logfile)