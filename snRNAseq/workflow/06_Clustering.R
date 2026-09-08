# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(ggplot2)

# Set seed for reproducibility
set.seed(100)

# Load snakemake variables
rds_input <- snakemake@input[["rds"]]
res <- snakemake@params[["resolution"]]
plot1 <- snakemake@output[["plot1"]]
plot2 <- snakemake@output[["plot2"]]
table_output <- snakemake@output[["table"]]
sample <- snakemake@params[["sample"]]

## Load data
sfi <- readRDS(rds_input)

## Print resolution
print(paste0(res, collapse = ", "))

# Assign identities to cells based on the present cluster
Idents(sfi) <- sfi[[]][, paste0("resolution_", res)]

# How many clusters do we have?
cli <- levels(sfi@active.ident) %>% length()

## Print UMAP
png(plot1)
DimPlot(sfi, label = TRUE) +
  ggtitle(sample) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Custom scoring metric (Roy Francis) which accounts for both adjusted p value AND log2FC
ranker <- function(m, pv = "p_val_adj", lf = "avg_log2FC") {
    m %>%
        mutate(score = -log10(!!sym(pv)) * abs(!!sym(lf))) %>%
        arrange(-score) %>%
        return()}

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

## Print dotplot
p <- DotPlot(sfi, features = unique(markers2plot$gene),
             cols = "RdBu") +
  coord_flip() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    axis.title = element_blank()
  )

ggsave(plot2, p, height = 3 + (as.integer(cli) * 0.6), width = 4 + (as.integer(cli) * 0.4), bg = "white")

# Save marker list
write_xlsx(x = markers, path = table_output)

# Close log file
sink(type = "message")
sink()
close(logfile)