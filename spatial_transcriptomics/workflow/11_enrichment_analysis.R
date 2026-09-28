###################### README
# Author: Inika Prasad, inika.prasad@gmail.com
# Date created: 2026-04-23
# Aim: To run GO term analysis on differentially expressed genes identified
# from the Border Zone (cluster 10) between vehicle and cCpE (after refinement)

# And visualize the results as heatmaps

# This is not a part of the Snakemake pipeline. Run independently

################### Load libraries #############################
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(tidyr)
library(stringr)
library(tibble)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)

set.seed(100)

##################### Filenames to save ##########################
logfile <- "../results/11_enrichment_analysis/BorderZone/BorderZone_logfile.xlsx"
gores_file <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOresults.xlsx"
goterm_groupfile <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOterm_groups.xlsx"
supp_file <- "../results/11_enrichment_analysis/BorderZone/Supplementary_BZ.xlsx"

heatmap_full <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmapfull.pdf"
heatmap_trimmed <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmaptrim.pdf"
heatmaptrimgenes_file <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmaptrimGENES.xlsx"

##################### Load and prep data ######################
df2 <- read_xlsx("../results/10_Cytoscape/BorderZone/noFCcutoff_intersections_ZANIA_ZAURAK.xlsx")
colnames(df2)

# Extract genenames & log2FC values
rankedlist <- df2$avg_log2FC
names(rankedlist) <- df2$gene_id

View(df2)

# Save DEGs for supplementary file
deg_df <- df2 %>%
  dplyr::select(c(gene, gene_id, p_val, avg_log2FC, pct.1, pct.2, p_val_adj))

##################### Run GO term analysis #####################
# unranked all
all_GO <- enrichGO(gene = names(rankedlist),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "ALL",
                   minGSSize = 3,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
dotplot(all_GO) + ggtitle("All DEGs GO ORA")

GOres <- all_GO@result
GOres %>% View()

# Save GO results
write_xlsx(GOres, gores_file)

################ Extract genes #################
# Extract each gene, and the GO terms its associated with
genes_long <- GOres %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest(geneID)

# How many genes are encompassed by GO terms? 49 genes
unique_genes <- unique(genes_long$geneID)
unique_genes

# Can we group the GO terms?
all_GO@result %>% dplyr::pull(Description) %>% View()
# Grouped GO Terms with help from Microsoft CoPilot, 2026-04-23
# And wrote the mapping in an excel file
# Load the mapping file:

mapfile <- read_xlsx(goterm_groupfile)
head(mapfile)

# Define the grouping
GOtermdf <- left_join(all_GO@result, mapfile, by = "Description")

# Make a supplementary file
write_xlsx(list(BorderZone_DEGs = deg_df,
                BorderZone_GO = GOtermdf),
           supp_file)

################# Map genes to Categories
View(genes_long)

nrow(genes_long)
nrow(mapfile) # 50 terms
intersect(genes_long$Description, mapfile$Description) # 50 terms

# Join info on genes & terms, with the categorization
gores_cat <- left_join(genes_long, mapfile, join_by("Description"))

# Join info on genes & terms & categorization with the DEGinfo
gores_cat <- left_join(gores_cat, deg_df, by = join_by(geneID == gene_id))

################ How many genes per category?
genes_per_category <- gores_cat %>%
  distinct(Curated_Category, geneID) %>%
  count(Curated_Category, name = "N_genes")

genes_per_category

################## How many genes in multiple categories?
genes_multiple_categories <- gores_cat %>%
  distinct(geneID, Curated_Category) %>%
  dplyr::count(geneID, name = "N_categories") %>%
  filter(N_categories > 1) %>%
  arrange(desc(N_categories))

genes_multiple_categories

#################### Order the categories for heatmaps
# This is just to order the categories!
cats2retain <- c("Muscle development & differentiation",
                 "Sarcomere organization",
                 "Cell adhesion",
                 "ECM & Collagen organisation",
                 "Cytoskeleton & Actin organization",
                 "Mitochondrial function & metabolism",
                 "Translation & Ribosome assembly",
                 "Other Functions")

gores_cat$Curated_Category <- factor(gores_cat$Curated_Category, levels = cats2retain)

#################### OPTION1: Make a full heatmap
currentplotname <- heatmap_full

################## OPTION2: Make a trimmed heatmap
cats2retain <- c("Muscle development & differentiation",
                 "Sarcomere organization",
                 "Cell adhesion",
                 "ECM & Collagen organisation",
                 "Cytoskeleton & Actin organization")

currentplotname <- heatmap_trimmed

################### Make the matrix for the heatmap
heatmap_matrix <- gores_cat %>%
  dplyr::filter(Curated_Category %in% cats2retain) %>%
  distinct(Curated_Category, geneID) %>%
  left_join(., deg_df, by = join_by(geneID == gene_id)) %>%
  dplyr::select(c(Curated_Category, geneID, avg_log2FC)) %>%
  pivot_wider(names_from = geneID,
              values_from = avg_log2FC,
              values_fill = list(avg_log2FC = NA)) %>%
  arrange(Curated_Category) %>%
  column_to_rownames("Curated_Category") %>%
  as.matrix()

heatmap_binary <- (heatmap_matrix != 0) * 1

first_one <- apply(heatmap_binary, 2, function(x) {
  i <- which(x == 1)
  if (length(i) == 0) Inf else i[1]
})

heatmap_matrix <- heatmap_matrix[, order(first_one)]

View(heatmap_matrix)

# Make the column names good
colnames(heatmap_matrix) <- str_to_sentence(colnames(heatmap_matrix))

################## Plot the heatmap #############
minval <- min(heatmap_matrix, na.rm = TRUE)
maxval <- max(heatmap_matrix, na.rm = TRUE)
minval
maxval
scale_end <- 2.5

if (scale_end > abs(minval) & scale_end > abs(maxval)) {
  pheatmap(
  heatmap_matrix,
  scale = "none",
  na_col = "white",
  breaks = seq(-scale_end, scale_end, length.out = 101),
  border_color = "lightgrey",
  filename = currentplotname,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  cellwidth = 12,
  cellheight = 25,
  width = 15
)
} else {
  print("pick better scale end! It needs to be bigger than the absolute values of
         the min and max of the heatmap")
}
