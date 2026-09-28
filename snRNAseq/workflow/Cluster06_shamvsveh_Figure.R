#################################### Load libraries
library(dplyr)
library(Seurat)
library(viridis)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(writexl)

set.seed(100)

#################################### Load Mapping File (Zeyu)
# To convert identifiers back to NCBI gene IDs
identifiers <- read.csv("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.tsv", sep = "\t")
identifiermap <- setNames(identifiers$gene_id, identifiers$genename_Hs_PwNameFormat)

###################################### Fig. 6
degfile <- "results/DGE/ALGORAB/ALGORABEPDC_A_DGE_results.xlsx"

my_palette <- colorRampPalette(plasma(20))(100) # custom viridis colors
options(enrichplot.colours = my_palette) # set viridis as colorscale

# Load DEG list
deglist <- read_xlsx(degfile, sheet = "sig_degs") %>%
    mutate(gene_id = identifiermap[gene]) # convert gene names to NCBI identifiers


# Extract genenames & log2FC values
filtlist <- setNames(deglist$avg_log2FC, deglist$gene_id)
length(filtlist) # 624 genes

# Divide list into up- and down-regulated genes
filtlist_up <- filtlist[filtlist > 0]
filtlist_down <- filtlist[filtlist < 0]

length(filtlist_up) #388 genes
length(filtlist_down) #236 genes

# unranked upregulated
upreg_GO <- enrichGO(gene = names(filtlist_up),
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     minGSSize = 3,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
as.data.frame(upreg_GO) %>% View()
dotplot(upreg_GO) + ggtitle("Upreg DEGs GO ORA")
# Save dotplot for upregulated


# unranked downregulated
downreg_GO <- enrichGO(gene = names(filtlist_down),
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       minGSSize = 3,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
as.data.frame(downreg_GO) %>% View()
dotplot(downreg_GO) + ggtitle("Downreg DEGs GO ORA")
# Save dotplot for downregulated

############################# Make supplementary file ############################
Supp_degdf <- deglist %>%
  dplyr::select(-c(score, score2, score3, score4, score5))

Supp_GOdf_up <- as.data.frame(upreg_GO) 

Supp_GOdf_down <- as.data.frame(downreg_GO) 

write_xlsx(list(Cluster06_DEGs_shamveh = Supp_degdf,
                Cluster06_downregGOinveh = Supp_GOdf_down,
                Cluster06_upregGOinveh = Supp_GOdf_up),
           path = "results/figures/Cluster06_supplementary.xlsx")

############################## Fig. 6d DOTPLOT EPDC-A sham vs. cryoinjury (veh) ###########

# Consider overlap between terms
intersect(upreg_GO@result$Description, downreg_GO@result$Description) # 22 overlaps

# Terms in upreg but not downreg
setdiff(upreg_GO@result$Description, downreg_GO@result$Description) %>% View()
terms2plot_up <- c("extracellular structure organization", # in both lists
                   "response to hypoxia", # no "hypoxia" in the other list
                   "transforming growth factor beta receptor signaling pathway", # no "transforming" in the other list
                   "regulation of anoikis", # no "anoikis" in the other list
                   "actin filament organization") # the other list has "actin filament polymerization (2 terms, while this list has 11 terms)"

# Terms in downreg but not upreg
setdiff(downreg_GO@result$Description, upreg_GO@result$Description) %>% View()
terms2plot_down <- c("extracellular structure organization", # in both lists
                     "collagen fibril organization", # no "collagen" in the other list
                     "regulation of cell junction assembly", #no "junction" in the other ist
                     "metallopeptidase activity", # no "metallopeptidase" in the other list
                     "Rho protein signal transduction") # no "rho" in the other list


my_palette <- colorRampPalette(viridis(20))(100) # custom viridis colors
options(enrichplot.colours = my_palette) # set viridis as colorscale

# SANITY CHECK: Plot terms of interest (raw dotplots)
p_up_og <- dotplot(upreg_GO, showCategory = terms2plot_up)
p_down_og <- dotplot(downreg_GO, showCategory = terms2plot_down)
p_up_og + p_down_og

# Extract results for custom dotplot
df_up <- upreg_GO@result %>%
  filter(Description %in% terms2plot_up) %>%
  mutate(Direction = "Upregulated")

df_down <- downreg_GO@result %>%
  filter(Description %in% terms2plot_down) %>%
  mutate(Direction = "Downregulated")

# Combine up- and down-reg terms results
df <- bind_rows(df_up, df_down) %>%
  # Parse GeneRatio (e.g. "10/203")
  tidyr::separate(GeneRatio, into = c("hits", "total"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio = hits / total) %>%
  dplyr::select(-hits, -total) %>%
  group_by(Direction) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  # mutate(Description = reorder_within(Description, GeneRatio, Direction)) %>%
  ungroup()

# Order by Gene Ratio
order <- df %>% arrange(GeneRatio) %>% dplyr::pull(Description) %>% unique()
df$Description <- factor(df$Description, levels = order)
df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))

# Make dotplot
ggplot(df,
       aes(x = GeneRatio,
           y = Description,
           colour = p.adjust,
           size = Count)) +
  geom_point() +
  facet_wrap(~Direction, ncol = 1, scales = "free_y") +
  scale_colour_gradientn(colours = my_palette, trans = "reverse") +
  theme_bw() +
  labs(title = "Cluster 06: Sham vs. Vehicle",
       x = "Gene Ratio",
       y = NULL,
       colour = "Adjusted p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "left",
    strip.background = element_rect(fill = "#1c368d"),
    strip.text = element_text(colour = "white", face = "bold"))

ggsave("results/figures/Cluster06_dotplot_shamveh.pdf", height = 3.5, width = 8)

################################ Make violin plots

################# Fig. 6e: Sham vs Veh: Collect genes for term of interest

# Find and specify genes of interest
toi <- "extracellular structure organization" # term of interest
goi <- df[df$Description == toi,]$geneID
genes <- strsplit(goi, "/") %>% unlist() # genes of interest

# Switch back to genename_Hs_Pw identifiers for plotting
genes_Pw <- identifiermap[identifiermap %in% genes] #names of this vector for plotting
genes2plot <- names(genes_Pw)
genes_Pw

{ # Epi peelout
input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
reduction <- "umap_peeled.harmony"
datasetname <- "10dpci_ctrl_Epipeelout"
clustering <- "merged_clusters"
clustering_infofile <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/10dpci_ctrl_EpiPeelout_CUSTOMres.xlsx"
df <- read_xlsx(clustering_infofile, sheet = clustering)
map <- setNames(df$Naming_manuscript, df$Cluster) # (desired naming, original naming)
cols <- setNames(df$cols, df$Naming_manuscript) # (colors, desired naming)
}

obj <- readRDS(input_rds)
obj2 <- JoinLayers(obj)
# make new metadata called named_clusters
obj2@meta.data$named_clusters <- factor(
  map[as.character(obj2@meta.data[[clustering]])],
  levels = unname(map))

clustering <- "named_clusters"
Idents(obj2) <- clustering
levels(obj2)

DimPlot(obj2, label = TRUE)

# Subset object to only include sham and vehicle
obj3 <- subset(obj2, subset = orig.ident %in% c("14dpsham", "10dpci_veh") & merged_clusters == "EPDC_A")

# Specify colors for the violin plot
cond_cols <- c("PM_naive" = "#b7fff9",
               "14dpsham" = "#ffacf8",
               "10dpci_veh" = "#b2ffab",
               "10dpci_cCpE" = "#fffbc0")

table(obj3$orig.ident, obj3$merged_clusters)


# Exploratory phase: Break up the list of top_genes into 20-gene chunks so printing is easy
gene_groups <- split(genes2plot, ceiling(seq_along(genes2plot) / 20))
pdf("results/figures/ECM_Cluster06_violinplots_upreg.pdf", height = 10, width = 10)
for (i in gene_groups) {
p <- VlnPlot(obj3,
             features = i,
             ncol = 5,
             alpha = 0,
             cols = cond_cols,
             group.by = "orig.ident")
  print(p)
}

dev.off()

genes2plot <- c("ADAMTS1", "ADAMTS4", "ADAMTS15", "MMP19-PW",
                     "COL1A1-PW", "COL3A1-PW", "COL6A1", "SERPINH1-PW")
# To use in the final figure
identifiermap[genes2plot]

# Final plots
VlnPlot(obj3,
        features = genes2plot,
        alpha = 0,
        group.by = "orig.ident",
        cols = cond_cols,
        ncol = 4) &
  scale_x_discrete(labels = c("Sham", "Vehicle")) &
  labs(x = NULL, y = NULL) &
  theme(legend.position = "none")
ggsave("results/figures/Cluster06_violins.pdf", width = 7, height = 4)