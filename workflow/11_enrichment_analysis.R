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

##################### Load and prep data ######################
df2 <- read_xlsx("../results/10_Cytoscape/BorderZone/noFCcutoff_intersections_ZANIA_ZAURAK.xlsx")

# Extract genenames & log2FC values
rankedlist <- df2$avg_log2FC
names(rankedlist) <- df2$genename_Hs

# sort in decreasing order for ranked analysis (not done here)
rankedlist <- sort(rankedlist, decreasing = TRUE)

##################### Filenames to save ##########################
logfile <- "../results/11_enrichment_analysis/BorderZone/BorderZone_logfile.xlsx"
gores_file <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOresults.xlsx"
goterm_groupfile <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOterm_groups.xlsx"

heatmap_full <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmapfull.pdf"
heatmap_trimmed <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmaptrim.pdf"
heatmaptrimgenes_file <- "../results/11_enrichment_analysis/BorderZone/BorderZone_GOtermgroups_heatmaptrimGENES.txt"
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
genes_long <- GOres %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  unnest(geneID)

# How many genes are encompassed by GO terms? 62 genes
unique_genes <- unique(genes_long$geneID)

# Can we group the GO terms?
all_GO@result %>% dplyr::pull(Description) %>% View()
# Grouped GO Terms with help from Microsoft CoPilot, 2026-04-23
# And wrote the mapping in an excel file
# Load the mapping file:

mapfile <- read_xlsx(goterm_groupfile)
head(mapfile)

################# Map genes to Categories
View(genes_long)

gores_cat <- left_join(genes_long, mapfile, join_by("Description"))
View(gores_cat)

################ How many genes per category?
genes_per_category <- gores_cat %>%
  distinct(Curated_Category, geneID) %>%
  count(Curated_Category, name = "N_genes")

genes_per_category

################## How many genes in multiple categories?
genes_multiple_categories <- gores_cat %>%
  distinct(geneID, Curated_Category) %>%
  count(geneID, name = "N_categories") %>%
  filter(N_categories > 1) %>%
  arrange(desc(N_categories))

genes_multiple_categories

#################### OPTION1: Make a full heatmap
# This is just to order the categories!
cats2retain <- c("Muscle_Development", "Sarcomere_Organization", "Heart_Contractility",
                 "Cell_Adhesion", "Collagen_Organization", "Actin_Cytoskeleton", "Mitochondrial",
                 "Translation_Ribosome", "Other_Functions")

gores_cat$Curated_Category <- factor(gores_cat$Curated_Category, levels = cats2retain)

heatmap_matrix <- gores_cat %>%
  distinct(Curated_Category, geneID) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = geneID,
              values_from = value,
              values_fill = list(value = 0)) %>%
  arrange(Curated_Category) %>%
  column_to_rownames("Curated_Category") %>%
  as.matrix()

View(heatmap_matrix)

# Make the column names good
colnames(heatmap_matrix) <- str_to_sentence(colnames(heatmap_matrix))

currentplotname <- heatmap_full

################## OPTION2: Make a trimmed heatmap
cats2retain <- c("Muscle_Development", "Sarcomere_Organization", "Heart_Contractility",
                 "Cell_Adhesion", "Collagen_Organization", "Actin_Cytoskeleton")

genes2exclude <- c("RPS10", "RPS16", "RPL22", "RPS19")

heatmap_matrix <- gores_cat %>%
  dplyr::filter(Curated_Category %in% cats2retain) %>%
  dplyr::filter(!geneID %in% genes2exclude) %>%
  distinct(Curated_Category, geneID) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = geneID,
              values_from = value,
              values_fill = list(value = 0)) %>%
  column_to_rownames("Curated_Category") %>%
  as.matrix()

# Save the gene names for use in Cytoscape
writeLines(colnames(heatmap_matrix), heatmaptrimgenes_file)

colnames(heatmap_matrix) <- str_to_sentence(colnames(heatmap_matrix))

View(heatmap_matrix)

################## Plot the heatmap #############
pheatmap(
  heatmap_matrix,
  scale = "none",
  color = c("white", "#9083C3"),
  border_color = "lightgrey",
  filename = currentplotname,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  cellwidth = 12,
  cellheight = 12,
  width = 15
)
