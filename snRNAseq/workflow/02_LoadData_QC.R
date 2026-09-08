###################################### Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

###################################### Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

#################################### Load Snakemake variables
raw_data <- snakemake@input[["raw_data"]]
sample <- snakemake@wildcards[["sample"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
plot3 <- snakemake@output[["plot3"]]
table_output <- snakemake@output[["table"]]
rds_output <- snakemake@output[["rds"]]

###################################### Load data
so <- Read10X_h5(raw_data)

###################################### Make Seurat object
su <- CreateSeuratObject(so, project = sample)

###################################### QC plot & save
VlnPlot(su, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, alpha = 0.1)
ggsave(filename = plot1, bg = "transparent")

###################################### Basic info
# Plot nFeature vs. nCount in cells
nCount_vs_nFeature <- FeatureScatter(su,
                                     feature1 = "nCount_RNA",
                                     feature2 = "nFeature_RNA") +
  ggtitle(sample)

ggsave(filename = plot2, bg = "transparent")


## How many cells express a given gene? Help decide min.cells threshold
su_countmat <- LayerData(su, assay = "RNA", layer = "counts")
gene_counts <- rowSums(su_countmat > 0)
png(plot3)
hist(gene_counts, breaks = 100, main = "Number of Cells Expressing Each Gene")
dev.off()

###################################### Compute QC metrics
su_summary <- su[[]] %>%
  summarise(
    cells = n(),
    total_features = dim(su)[1],
    mean_nCount = mean(nCount_RNA, na.rm = TRUE),
    mean_nFeatures = mean(nFeature_RNA, na.rm = TRUE)
  )

###################################### Save QC metrics
write_xlsx(x = su_summary, path = table_output)

###################################### Save RDS object
saveRDS(object = su, file = rds_output)

###################################### Close log file
sink(type = "message")
sink()
close(logfile)