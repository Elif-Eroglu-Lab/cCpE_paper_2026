########################################## README
# Author: Inika Prasad
# Date: June 2026
# Aim: to identify genes which follow a pattern of
# highly expressed in 10,
# expressed in 02 and 09 upon injury

########################################### Load libraries
library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(writexl)
library(tibble)
library(readxl)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(viridis)
library(limma)
library(stringr)
library(purrr)

set.seed(100)

source("workflow/functions.R")

########################################### Read in Seurat Object

{ #10dpci_ctrl Epi Peelout: get the dataset
input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
reduction <- "umap_peeled.harmony"
celltypes <- "merged_clusters"

obj <- readRDS(input_rds)
Idents(obj) <- celltypes
obj2 <- JoinLayers(obj)

obj2$Condition_Cluster <- paste0(obj2$orig.ident, "_", obj2$merged_clusters)

# Get the colors
clustering_infofile <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/10dpci_ctrl_EpiPeelout_CUSTOMres.xlsx"
clustering <- "merged_clusters"
df <- read_xlsx(clustering_infofile, sheet = clustering)
cols <- df$cols
names(cols) <- df$Cluster

# Set levels
obj2$merged_clusters <- factor(obj2$merged_clusters, df$Cluster)
Idents(obj) <- celltypes
}

################### Define filenames ############################################
mybasename <- "results/regen_gene_module/13_"

seedgenes_image <- paste0(mybasename, "seedgenes_module.png")
poi_heatmap <- paste0(mybasename, "patternofinterest.pdf")
poi_heatmap_scaled <- paste0(mybasename, "patternofinterest_scale.pdf")

corrhist <- paste0(mybasename, "corr_histogram.png")
corrhist2 <- paste0(mybasename, "corr_histogram2.pdf")
corrdf2save <- paste0(mybasename, "corr_seed_genes.xlsx")

modulescore_umap <- paste0(mybasename, "top_genes_MODULE_split.pdf")

topgenes_human_file <- paste0(mybasename, "top_genes_human_orthologues.txt")
topgenes_human_excel <- paste0(mybasename, "top_genes_human_orthologues.xlsx")

topgenes_GOexcel <- paste0(mybasename, "top_genes_GO.xlsx")
basicGOdotplot <- paste0(mybasename, "basicGOdotplot.png")
GOobject <- paste0(mybasename, "top_genes_GO.rds")

# Annotated GO Terms file
topgenes_GOexcel_ann <- paste0(mybasename, "top_genes_GO_annotated_by_inika.xlsx")
customGOdotplot <- paste0(mybasename, "customGOdotplot.pdf")
cnetplot_file <- paste0(mybasename, "GOterms_cnetplot.pdf")
cnetplot_edges <- paste0(mybasename, "cnetplot_edges.csv")
cnetplot_nodes <- paste0(mybasename, "cnetplot_nodes.csv")

aggexp_heatmap <- paste0(mybasename, "topgenesHEATMAP.png")

ccpe_heatmap_mesenchclusts <- paste0(mybasename, "topgenesHEATMAP_wccpe_mesenchclusts.pdf")
ccpe_heatmap_mesenchclusts_long <- paste0(mybasename, "topgenesHEATMAP_wccpe_mesenchclustslong.pdf")

ccpe_heatmap <- paste0(mybasename, "topgenesHEATMAP_wccpe.pdf")

cameraPR_res <- paste0(mybasename, "cameraPR_results.txt")

full_heatmap <- paste0(mybasename, "topgenes_HEATMAPallclusts.pdf")
almost_full_heatmap <- paste0(mybasename, "topgenes_HEATMAPalmostallclusts.pdf")

heatmap_wccpe <- paste0(mybasename, "min_interestingpattern_wccpespb.pdf")
heatmap_wccpe_dense <- paste0(mybasename, "min_interestingpattern_wccpespb_dense.pdf")
med_heatmap <- paste0(mybasename, "min_interestingpattern_someclustspb.pdf")
epdcs_heatmap <- paste0(mybasename, "min_interestingpattern_epdcclustspb.pdf")

med_heatmap_sorted <- paste0(mybasename, "min_interestingpattern_someclustspb_sorted.pdf")
med_heatmap_long_sorted <- paste0(mybasename, "min_interestingpattern_someclustspb_longsorted.pdf")

supp_file <- paste0(mybasename, "regen_module_supplementary_genes.xlsx")

############################ Remove groups with too few cells ####################
{
# Identify cluster & condition groups with <3 cells
freq_tab <- table(obj2$orig.ident, obj2$merged_clusters) %>%
  as.data.frame() %>%
  rename(., Condition = Var1) %>%
  rename(., Cluster = Var2) %>%
  mutate(Condition_Cluster = paste0(Condition, "_", Cluster))

torem <- freq_tab %>% dplyr::filter(Freq < 3)
torem

# Remove groups with < 3 cells
seu_data <- subset(obj2,
                   subset = !(obj2$Condition_Cluster %in% torem$Condition_Cluster))


# Check subsetting was proper
seu_data$orig.ident %>% unique()
table(seu_data$orig.ident, seu_data$merged_clusters)

}

# Plots to check subsetting was proper
DimPlot(seu_data, label = TRUE, pt.size = 2.5, group.by = "orig.ident")
DimPlot(seu_data, label = TRUE, pt.size = 2.5)
dim(seu_data)

##################### Pseudobulk data ################################
pb <- AggregateExpression(seu_data,
                          assays = "RNA",
                          group.by = c("orig.ident", "merged_clusters"),
                          return.seurat = TRUE)

##################### Investigate expression of key genes ############

# Read markers from excel file
markers <- read_xlsx("results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/custom_res1_markers.xlsx")
markers$cluster <- factor(markers$cluster, levels = names(cols))
markers$cluster %>% unique()

# Extract Cluster 10 (inj-EPDCs) top markers
genes2plot <- markers %>%
  dplyr::filter(cluster == 10) %>%
  arrange(cluster, -(abs(avg_log2FC))) %>%
  slice_head(n = 50) %>%
  pull(gene)

genes2plot

##################### Make metadata for pseudobulked groups ############
# Extract expression values
expr <- pb@assays$RNA$scale.data

# Wrange metadata into shape
my_metadata <- pb@meta.data %>% rownames_to_column("group")
my_metadata <- my_metadata %>%
  mutate(orig.ident = gsub("\\_.*", "", orig.ident)) %>%
  mutate(orig.ident = gsub("g", "", orig.ident)) %>%
  mutate(merged_clusters = gsub("-", "_", merged_clusters))

# Levels of metadata
my_metadata$merged_clusters <- factor(my_metadata$merged_clusters, levels = names(cols))
my_metadata$orig.ident <- factor(my_metadata$orig.ident, levels = c("PM-naive", "14dpsham", "10dpci-veh", "10dpci-cCpE"))

# Arrange by levels of metadata
my_metadata <- my_metadata %>%
  arrange(orig.ident, merged_clusters)

# ############ Define gene pattern of interest with average expression ##################################
# Subset the dataset to remove the groups not needed to define the pattern of interest
g2retain <- my_metadata %>%
  dplyr::filter(orig.ident != "10dpci-cCpE") %>%
  dplyr::filter(merged_clusters %in% c("02", "09", "10")) %>%
  dplyr::pull(group)

g2retain

expr_sub <- expr[, g2retain]

# # Define the expression profile of the seed genes
seed_genes <- c("TNC-PW", "CHST15", "TGM2", "OXTR", "TGFBI", "CILP",
                "SCN5A", "SLC2A1-PW", "SLC1A5", "CYSLTR1-PW", "CPM-PW")

# seed_genes <- c("OXTR", "CILP", "SLC2A1-PW")

# seed_genes <- genes2plot

seed_profile <- colMeans(expr_sub[seed_genes, , drop = FALSE])
seed_profile %>% View()

############## Plot seed profile and genes for pattern ########################
# Make the column annotation table
mat_annot <- my_metadata %>%
  column_to_rownames("group") %>%
  mutate(Condition = orig.ident,
         Cluster = merged_clusters) %>%
  dplyr::select(c(Condition, Cluster))

# Define color palette
my_palette <- colorRampPalette(plasma(20))(100)

# Specify colors for the column annotations
ann_colors <- list(Condition = c("PM-naive" = "#b7fff9",
                                 "14dpsham" = "#ffacf8",
                                 "10dpci-veh" = "#b2ffab",
                                 "10dpci-cCpE" = "#fffbc0",
                                 "14dpci" = "#cbabff"),
                   Cluster = cols)

# Compile pattern of interest and component genes
Mgenes <- expr_sub[seed_genes, g2retain]
Pattern_of_interest <- seed_profile[g2retain]

M <- rbind(Mgenes, Pattern_of_interest)
rownames(M) <- str_to_sentence(rownames(M))

# Specify gaps in the columns
colgaps <- c(2, 5)
rowgaps <- length(seed_genes)

# Make the heatmap with seed genes and pattern of interest
pheatmap::pheatmap(
  M,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  border_color = NA,
  color = my_palette,
  annotation_col = mat_annot,
  annotation_colors = ann_colors,
  # legend_breaks = c(-1, 0, 1),
  gaps_col = colgaps,
  gaps_row = rowgaps,
  filename = poi_heatmap,
  scale = "none",
  # filename = poi_heatmap_scaled,
  # scale = "row",
  show_rownames = TRUE,
  cellheight = 12,
  cellwidth = 14,
  width = 8,
  angle_col = 45)

dev.off()


############ Correlate genes with seed gene profile ###########################
# Correlate genes with the seed gene profile
gene_corr <- apply(expr_sub,
                   1, # for every row (i.e. gene)
                   function(x) cor(x, seed_profile, method = "spearman"))

# Rank genes by highest correlation
corr_df <- data.frame(gene = names(gene_corr),
                      correlation = gene_corr)

corr_df <- corr_df[order(-corr_df$correlation), ]

# Save excel file of correlation coefficients
write_xlsx(corr_df, corrdf2save)

# Save histogram of correlation coefficients
png(corrhist)
hist(corr_df$correlation)
dev.off()

View(corr_df)

############# Fetch top correlated genes
corr_df <- read_xlsx(corrdf2save)

# Check how many NA values we have
is.na(corr_df$correlation) %>% sum() #4053 genes have NA values

# Make a publication quality histogram
ggplot(corr_df, aes(x = correlation)) +
  geom_histogram(color = "white", fill = "#868686",
                 binwidth = 0.1, boundary = 0, closed = "right") +
  ylab("No. of genes") +
  xlab("Spearman's correlation coefficient") +
  theme_classic()

ggsave(corrhist2, height = 4, width = 4)

hist(corr_df$correlation)
corr_df$correlation_og <- corr_df$correlation
corr_df$correlation[is.na(corr_df$correlation)] <- 0
top_genes <- corr_df$gene[corr_df$correlation > 0.9]
top_genes # 149 genes

############### Module expression of these genes

obj2 <- AddModuleScore(object = obj2, features = list(top_genes), name = "top_genes")
obj2$top_genes1

p <- FeaturePlot(obj2,
            pt.size = 2,
            reduction = reduction,
            features = "top_genes1",
            split.by = "orig.ident",
            keep.scale = "feature",
            cols = c("#cbf9ff", "#0400ff"),
            order = TRUE) &
    theme(legend.position = "right")
p

ggsave(modulescore_umap, width = 38)


############### GO enrichment of these genes
# get human names for the GO terms
# Load Zeyu's file
identifiers <- read.csv("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.tsv", sep = "\t")
colnames(identifiers)

identifiers <- identifiers %>%
  dplyr::select(c("genename_Hs_PwNameFormat", "genename_Hs"))

# Get human names for top genes
topgenes_human_df <- identifiers %>%
  dplyr::filter(genename_Hs_PwNameFormat %in% top_genes)
  
topgenes_human <- topgenes_human_df %>%
  dplyr::pull(genename_Hs) %>%
  unique()

# Run GO term enrichment analysis
all_GO <- enrichGO(gene = topgenes_human,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "ALL",
                   minGSSize = 3,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
as.data.frame(all_GO) %>% View()
dotplot(all_GO) + ggtitle("Regen module GO terms")
ggsave(basicGOdotplot)

# Save human gene names for top genes
writeLines(topgenes_human, topgenes_human_file)
write_xlsx(topgenes_human_df, topgenes_human_excel)

# Save GO results excel file
write_xlsx(as.data.frame(all_GO), topgenes_GOexcel)

# Save GO results analysis
saveRDS(all_GO, GOobject)
# all_GO <- readRDS(GOobject)

#################### Make custom dotplot
my_palette <- colorRampPalette(viridis(20))(100) %>% rev() # custom viridis colors
options(enrichplot.colours = my_palette) # set viridis as colorscale

all_GO <- readRDS(GOobject)
all_GOdf <- read_xlsx(topgenes_GOexcel_ann) %>% dplyr::filter(Keep == "Yes")
all_GOdf
dotplot(all_GO, showCategory = all_GOdf$Description, label_format = 40)
ggsave(customGOdotplot, height = 4, width = 7)

nrow(all_GO@result)

#################### Make cnetplot and export for cytoscape
# Make basic cnetplot to get a vibe check
p <- cnetplot(all_GO,
         showCategory = nrow(all_GO@result),
         color_item = "#B8B8FF",
         color_category = "#AA3E98",
         color_edge = "#797979")
p
ggsave(cnetplot_file, p, height = 15, width = 15)

# Save edge info
edges <- all_GO@result %>%
  # select(ID, Description, geneID, Count, p.adjust, ) %>%
  separate_rows(geneID, sep = "/") %>%
  rename(source = Description,
         target = geneID) %>%
  mutate(target_display = str_to_sentence(target),
         log10adj_pval = -log10(p.adjust),
         termorgene = log10adj_pval/log10adj_pval)
write.csv(edges, cnetplot_edges, row.names = FALSE, quote = FALSE)

# Save node info
terms <- data.frame(id = unique(edges$source), type = "Term")
genes <- data.frame(id = unique(edges$target), type = "Gene")
nodes <- bind_rows(terms, genes)
write.csv(nodes, cnetplot_nodes, row.names = FALSE, quote = FALSE)

############### Make a heatmap of pseudobulked data
filename <- aggexp_heatmap

pheatmap::pheatmap(
  expr[top_genes, ],
  # cluster_cols = FALSE,
  cluster_rows = FALSE,
  filename = filename,
  cellheight = 12,
  cellwidth = 14,
  angle_col = 45,
  width = 15,
  scale = "none")

dev.off()

#############################################################################
################# Plot top genes in the wider datasets
 #10dpci_ctrl Epi Peelout: get the dataset with ccpe

{ # see if these genes follow the pattern
  g2plot <- g2retain
  cellheight <- 1.2
  filename2 <- paste0(mybasename, "top_genes_coi.pdf")
}

{ # OPTION: Heatmap with all conditions, all clusters
  g2plot <- mat_annot %>% rownames()
  filename2 <- full_heatmap
  cellheight <- 1.2
}

{ # OPTION: Heatmap with all conditions and 06, 07, 08, 09. Squished
  g2plot <- mat_annot %>% dplyr::filter(Cluster %in% c("02", "09", "10", "EPDC_A")) %>% rownames()
  filename2 <- ccpe_heatmap_mesenchclusts
  cellheight <- 12
  gaps_cols <- c(3, 7, 11)
}

{ # OPTION: Heatmap with all conditions and 06, 07, 08, 09. Squished
  g2plot <- mat_annot %>%
    arrange(Cluster, Condition) %>%
    dplyr::filter(Cluster %in% c("02", "09", "10", "EPDC_A")) %>% rownames()
  filename2 <- ccpe_heatmap_mesenchclusts
  cellheight <- 1.2
  gaps_cols <- c(4, 8, 11)
}

############## Make the heatmap
my_palette <- colorRampPalette(plasma(20))(100)

M <- expr[top_genes, g2plot] # Extract top genes

# expr_scale <- pb@assays$RNA$data
# M <- expr[top_genes, g2plot] # Extract top genes

# expr_counts <- pb@assays$RNA$counts
# M <- expr_counts[top_genes, g2plot] # Extract top genes

# filename2 <- paste0(mybasename, "0.9topgenes_HEATMAPalmostallclusts_UNSCALED_ROWCLUST.pdf")

# Make the heatmap
pheatmap::pheatmap(
  M,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  color = my_palette,
  annotation_col = mat_annot,
  annotation_colors = ann_colors,
  # gaps_col = gaps_cols,
  filename = filename2,
  # treeheight_row = 200,
  show_rownames = FALSE,
  scale = "none",
  cellheight = cellheight,
  cellwidth = 14,
  width = 15,
  angle_col = 45)

dev.off()

######## Go make heatmap in figures.R
genes2plot <- top_genes

####################### Make comparisons between veh & cCpE
# Create metadata
Idents(obj2) <- "Condition_Cluster"

# Make sure data looks reasonable with condition_cluster
DimPlot(obj2)
Idents(obj2) %>% unique()

######## See where the genes end up in the DEGs

# Make the comparison sets (same cluster, across conditions)
df_wide <- mat_annot %>%
  mutate(Condition = gsub("-", "_", Condition)) %>%
  dplyr::select(Cluster, Condition) %>%
  mutate(Condition_Cluster = paste0(Condition, "_", Cluster)) %>%
  pivot_wider(
    names_from = Condition,
    values_from = Condition_Cluster
  )

# Make comparisons dataframe
comps4barcode <- data.frame(ident.2 = df_wide$`10dpci_veh`, ident.1 = df_wide$`10dpci_cCpE`) # make columns with gorups 2 compare
comps4barcode$cluster <- gsub("^(14dpsham_|10dpci_veh_)", "", comps4barcode$ident.2) # Create a column with Cluster info

# Make labels for comparisons
comps4barcode <- comps4barcode %>%
  mutate(comp = case_when(
         grepl("sham", ident.2) & grepl("veh", ident.1)   ~ "sham_veh",
         grepl("veh", ident.2) & grepl("cCpE", ident.1) ~ "veh_cCpE",
         TRUE ~ NA_character_)) %>%
  mutate(label = paste0("Cluster", cluster, "_", comp)) %>%
  drop_na()

comps4barcode

sink(cameraPR_res)
print(comps4barcode)

# MAY NEED TO RESET WORM.MAX
for (i in 1:nrow(comps4barcode)) {
  # print(i)
  print(comps4barcode[i,]$label)

  FCdf <- FoldChange(obj2,
                     ident.1 = comps4barcode[i,]$ident.1,
                     ident.2 = comps4barcode[i,]$ident.2)
  label <- comps4barcode[i,]$label

  # Arrange log2FC values in decreasing order
  FCdf <- FCdf %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column("gene")

  # Create a vector with log2FC values and genes
  ranks <- FCdf$avg_log2FC
  names(ranks) <- FCdf$gene

  # Find indices of genes in top_genes
  geneSetIndex <- which(names(ranks) %in% top_genes)
  geneSetIndex %>% length()

  # add info about which indices are in top_genes into dataframe
  FCdf$top_gene <- FCdf$gene %in% names(ranks[geneSetIndex])
  FCdf <- FCdf %>% rownames_to_column("top_gene_rank")

  # Add info with human gene names
  FCdf <- left_join(FCdf, identifiers, join_by(gene == genename_Hs_PwNameFormat))

  # View(FCdf)

  # Save data on fold changes, and which genes were in the top genes
  write_xlsx(list(topgenes = FCdf %>% dplyr::filter(top_gene == TRUE),
                  allgenes = FCdf),
            paste0(mybasename, "barcodeplots/barcodes_", label, ".xlsx"))

  barcodepng <- paste0(mybasename, "barcodeplots/", label, ".pdf")

  # Save the barcode plot
  pdf(barcodepng, height = 7, width = 7)
  barcodeplot_mod(
    ranks,
    quantile = c(-0, 0),
    worm = TRUE,
    worm.max = 3.5,
    rectangle1col = "#5891db",
    rectangle2col = "#ee3e35",
    xlab = "Log2 Fold Change",
    labels = c("Downregulated", "Upregulated"),
    index = geneSetIndex)
  dev.off()



  # Calculate if the genes of the set are highly ranked compared
  # to the background genes
  print(cameraPR(FCdf$avg_log2FC, geneSetIndex))
}

sink()

################################# Check for overlap with activated Epi DEGs
{
dge_df <- read_xlsx("results/DGE/ALGOL/ALGOL04_DGE_results.xlsx", sheet = "sig_degs")
ranked_dge_df <- dge_df %>%
  arrange(-avg_log2FC) %>%
  rownames_to_column("rank")

label <- "ALGOL04"

intersect(top_genes, dge_df$gene)

}

################################ Identify which genes are ERK targets
ERK_targets <- read.csv("results/testerdump/A Compendium of ERK targets.csv")

View(ERK_targets)

ERK_targets %>%
  dplyr::filter(Gene.Name %in% topgenes_human) %>% #Many rows per gene because different binding sites
  dplyr::pull(Gene.Name) %>%
  unique() %>%
  str_to_sentence()


#################################### Make supplementary file

comps <- c("Cluster 06" = "ClusterEPDC_A_veh_cCpE",
           "Cluster 07" = "Cluster09_veh_cCpE",
           "Cluster 08" = "Cluster10_veh_cCpE",
           "Cluster 09" = "Cluster02_veh_cCpE")

dflist <- list()

for (i in 1:length(comps)) {
  comp <- comps[i]
  compname <- names(comps)[i]
  compfile <- paste0(mybasename, "barcodeplots/barcodes_", comp, ".xlsx")

  file.exists(compfile) %>% print()

  # read in the file and get the columns you need
  compdf <- read_xlsx(compfile, sheet = "allgenes") %>%
    dplyr::select(c(gene, top_gene, genename_Hs, top_gene_rank, avg_log2FC))

  colnames(compdf)[colnames(compdf) == "top_gene_rank"] <- paste0(compname, ": rank")
  colnames(compdf)[colnames(compdf) == "avg_log2FC"] <- paste0(compname, ": avg log2FC")
  colnames(compdf)[colnames(compdf) == "top_gene"] <- "regenerative gene module"

  print(head(compdf))
  dflist[[i]] <- compdf
  }

dflist

library(purrr)
clustDFs <- dflist %>% purrr::reduce(full_join, by = c("gene", "regenerative gene module", "genename_Hs"))
dim(clustDFs)
View(clustDFs)

corrdf <- read_xlsx(corrdf2save)
dim(corrdf)

supp_table <- full_join(corrdf, clustDFs, by = "gene")
View(supp_table)

write_xlsx(supp_table, supp_file)

################### How do the original seed genes look?

supp_table %>% dplyr::filter(gene %in% seed_genes) %>% View()

#################### Look at fold changes for a random list vs. this list
# Choose a list of random 243 genes
rownames(expr)
corr_df


randgenes <- sample(corr_df$gene, size = length(top_genes), replace = FALSE)


clust09 <- FoldChange(
  obj2,
  ident.1 = "10dpci_veh_09",
  ident.2 = "10dpci_cCpE_09",
  # features = randgenes
 features = top_genes
)

clust02 <- FoldChange(
  obj2,
  ident.1 = "10dpci_veh_02",
  ident.2 = "10dpci_cCpE_02",
  # features = randgenes
 features = top_genes
)

clust10 <- FoldChange(
  obj2,
  ident.1 = "10dpci_veh_10",
  ident.2 = "10dpci_cCpE_10",
  # features = randgenes
 features = top_genes
)

x1 <- clust02
x2 <- clust10

plot(x1$avg_log2FC, x2$avg_log2FC)
abline(a = 0, b = 1, col = "red", lwd = 2)

# B1 = 0.26, Clust 02 has attenuated response relative to Clust 09
fit <- lm(x1$avg_log2FC ~ x2$avg_log2FC)
summary(fit)
abline(fit, col = "blue", lwd = 2)
t.test(x1$avg_log2FC, x2$avg_log2FC, paired = TRUE)

# x has larger mean than y
t.test(x = abs(x1$avg_log2FC), y = abs(x2$avg_log2FC), paired = TRUE, alternative = "greater")
t.test(x = x1$avg_log2FC, y = x2$avg_log2FC, paired = TRUE, alternative = "greater")


###########################
# Do the DGE testing
m <- FindMarkers(obj2,
                ident.1 = "10dpci_cCpE_09", #tst
                ident.2 = "10dpci_veh_09", #ref
                features = top_genes,
                group.by = "Condition_Cluster",
                test.use = "MAST",
                only.pos = FALSE,
                # min.pct = 0.25,
                logfc.threshold = 0,
                min.cells.group = 3,
                random.seed = 100,
                assay = "RNA",
                slot = "data",
                latent.vars = c("nFeature_RNA", "nCount_RNA")) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, everything()) %>%
  ranker()

View(m)

####################################### Comparison with old data

# genedf <- read_xlsx("results/regen_gene_module/Attempt09_dataslotSeurat/09_corr_seed_genes.xlsx")
# genedf <- read_xlsx("results/regen_gene_module/Attempt08_log1p/08_corr_seed_genes.xlsx")
# genedf <- read_xlsx("results/regen_gene_module/10_corr_seed_genes.xlsx")
genedf <- read_xlsx("results/regen_gene_module/Attempt11_dataslotSeurat_11genes/11_corr_seed_genes.xlsx")


genedf$correlation_og <- genedf$correlation
genedf$correlation[is.na(genedf$correlation)] <- 0
top_genes_old <- genedf$gene[genedf$correlation > 0.9]
top_genes_old

intersect(top_genes, top_genes_old)

#############################################

# Make function to run cameraPR many times on random & selected subset
run_camera_null <- function(ranks,
                            top_genes,
                            n_iter = 1000,
                            subset_size = 100,
                            seed = 1) {

set.seed(seed)
all_genes <- names(ranks)

# CHECKPOINT: Only compare against top_genes that are actually present in the ranked universe
top_genes_present <- intersect(top_genes, all_genes)

# CHECKPOINT: check subset size is smaller than gene set
if (length(top_genes_present) < subset_size) {
  stop("top_genes has fewer genes than subset_size after intersecting with ranks")
}

# Create vectors to populate with p-values
random_pvals  <- numeric(n_iter)
topgene_pvals <- numeric(n_iter)

for (i in seq_len(n_iter)) {

  # Negative control: fully random gene subset from the whole background
  rand_genes <- sample(all_genes, subset_size)
  rand_index <- which(all_genes %in% rand_genes)
  rand_res   <- cameraPR(ranks, rand_index)
  random_pvals[i] <- rand_res$PValue

  # Random subset drawn only from your top_genes list
  tg_genes <- sample(top_genes_present, subset_size)
  tg_index <- which(all_genes %in% tg_genes)
  tg_res   <- cameraPR(ranks, tg_index)
  topgene_pvals[i] <- tg_res$PValue
  }

  # Output dataframe
  data.frame(iteration    = seq_len(n_iter),
             random_pval  = random_pvals,
             topgene_pval = topgene_pvals)
}


################################
Idents(obj2) <- "Condition_Cluster"



################################
# Run for many comparisons
for (i in 1:nrow(comps4barcode)) {
  # print(i)
  print(comps4barcode[i,]$label)

  FCdf <- FoldChange(obj2,
                     ident.1 = comps4barcode[i,]$ident.1,
                     ident.2 = comps4barcode[i,]$ident.2)
  label <- comps4barcode[i,]$label

  # Arrange log2FC values in decreasing order
  FCdf <- FCdf %>%
    arrange(-avg_log2FC) %>%
    rownames_to_column("gene")

  # Create a vector with log2FC values and genes
  ranks <- FCdf$avg_log2FC
  names(ranks) <- FCdf$gene

# Run simulation of random genes vs. top genes
df_vals <- run_camera_null(ranks,
                top_genes,
                n_iter = 1000,
                subset_size = 100,
                seed = 1)

View(df_vals)

# -log10 of p-values to make plotting and testing easier
df_vals <- df_vals %>%
  mutate(random_log  = -log10(random_pval),
         topgene_log = -log10(topgene_pval))


# Make dataframe good for plotting
plot_df <- df_vals %>%
  dplyr::select(iteration, random_log, topgene_log) %>%
  tidyr::pivot_longer(cols = c(random_log, topgene_log),
                      names_to = "group",
                      values_to = "neglog10p") %>%
  mutate(group = recode(group,
                         random_log  = "100 Random genes",
                         topgene_log = "100 Module genes"))

# Make density plot
ggplot(plot_df, aes(x = neglog10p, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#FFB8D1", "#1446A0")) +
  labs(x = "-Log10(cameraPR p-value)",
       y = "Density",
       title = label) +
  theme_classic()

ggsave(paste0(mybasename, label, "randomgene_comparison.pdf"), height = 5, width = 6)
}

# Is testing necessary or worth it?
# Non-parametric comparison of distributions (robust to skew, doesn't assume normality)
wilcox.test(df_vals$topgene_log, df_vals$random_log)

# Optional: test whether the whole shape differs, not just central tendency
ks.test(df_vals$topgene_log, df_vals$random_log)
