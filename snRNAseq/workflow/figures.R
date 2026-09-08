################################ README
# Figures for veh-ccpe manuscript


################################ Load Libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(writexl)
library(readxl)
library(patchwork)
library(viridis)
library(pals)
library(scales)
library(pheatmap)
library(clusterProfiler)


set.seed(100)

source("workflow/functions.R")

################################ Transparency
# This will be overridden if you use another theme in your plot
theme_set(theme_classic() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA)
    )
)

################################ Select data
{ # 10dpci_ctrl
input_rds <- "results/12_integrated_DimReduction/10dpci_ctrl/DRscrub/10dpci_ctrl_DRscrub_dimred.rds"
datasetname <- "10dpci_ctrl"
reduction <- "umap.harmony"
clustering <- "harmony_resolution_0.6"
clustering_infofile <- "results/testerdump/10dpci_ctrl_harmony_resolution_0.6.xlsx"

ptsize <- 0.6

df <- read_xlsx(clustering_infofile)

map <- df$Naming_manuscript
names(map) <- df$Cluster

cols <- df$cols
names(cols) <- df$Naming_manuscript

}

{ # Epi peelout
input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
reduction <- "umap_peeled.harmony"
datasetname <- "10dpci_ctrl_Epipeelout"
clustering <- "merged_clusters"
# clustering <- "harmony_peeled_resolution_1.4"
ptsize <- 2
clustering_infofile <- "results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/10dpci_ctrl_EpiPeelout_CUSTOMres.xlsx"

df <- read_xlsx(clustering_infofile, sheet = clustering)

map <- df$Cluster
names(map) <- df$Cluster

cols <- df$cols
names(cols) <- df$Cluster

cols %>% show_col()
}

############################## Load data
obj <- readRDS(input_rds)
obj2 <- JoinLayers(obj)

# make new metadata called named_clusters
obj2@meta.data$named_clusters <- factor(
  map[as.character(obj2@meta.data[[clustering]])],
  levels = unname(map))
clustering <- "named_clusters"

Idents(obj2) <- clustering
levels(obj2)

show_col(cols)

################################ Make UMAP

DimPlot(obj2, reduction = reduction,
        group.by = clustering,
        label = TRUE,
        cols = cols,
        # label.size = ptsize,
        repel = TRUE,
        pt.size = ptsize) +
  labs(x = "UMAP1", y = "UMAP2") &
  theme(legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent",
                                   color = NA)
  )

ggsave(paste0("results/figures/", datasetname, "_namedUMAP.png"), width = 10)

ggsave(paste0("results/figures/", datasetname, "_EpiPeelout_pre-merging.png"), width = 15, height = 12)
ggsave(paste0("results/figures/", datasetname, "_bigUMAPlabelled.png"), width = 12)


##################### Get color legend
# Convert to data frame
df <- data.frame(
  name = names(cols),
  color = unname(cols),
  x = seq_along(cols),
  y = 1
)

ggplot(df, aes(x = -x, y = y, fill = color)) +
  geom_tile(width = 0.9, height = 0.3) +  # skinny rectangles
  scale_fill_identity() +                 # use actual color values
  geom_text(aes(label = name), hjust = -0.8, vjust = 0, size = 12) +
  theme_void() +
  theme(legend.position = "none") +
  coord_flip()

filename <- "results/figures/clusterkey_10dpci_ctrl.pdf"
ggsave(filename, height = nrow(df), width = 2)

################################ Split by condition
# Extended Data Fig. 5a, Extended Data Fig. 6b

table(obj2$orig.ident)

plotlist <- list()

for (cond in unique(obj2$orig.ident)) {
  # make a new column with cluster identities same as usual
  obj2@meta.data[["condition_clust"]] <- obj2@meta.data[[clustering]]

  # Add "Other" as a level into metadata so you don't get NAs
  levels(obj2@meta.data$condition_clust) <- c(
    levels(obj2@meta.data$condition_clust),
    "Other"
  )

  # for every cell not in the condition of interest,
  # make the new column "Other"
  obj2@meta.data[["condition_clust"]][obj2$orig.ident != cond] <- "Other"

  p <- DimPlot(obj2, reduction = reduction,
          # label = TRUE,
          cols = c(cols, "Other" = "#c5c5c5"),
          order = c(names(cols), "Other"),
          group.by = "condition_clust",
          label.size = 6,
          repel = TRUE,
          pt.size = ptsize) +
    labs(x = "UMAP1", y = "UMAP2") +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(cond)

  plotlist[[cond]] <- p

}

(plotlist[[4]] + plotlist[[3]] + plotlist[[1]] + plotlist[[2]]) + plot_layout(nrow = 1)
(plotlist[[1]] + plotlist[[2]] + plotlist[[3]] + plotlist[[4]]) + plot_layout(nrow = 1)
ggsave(paste0("results/figures/", datasetname, "bycondition.png"), width = 30, height = 8)


################################ Specify genes for feature plots
{ # Post-Epi peelout: Fig. 5e
  genes2plot <- c("WT1", "TBX18", "TCF21", "BNC1",
                  "LRRN4", "PODXL", "CLDN15", "CLDN6-LOC138304364", "CLDN7",
                  "TNNT2", "TNNI3", "PECAM1", "EGFL7",
                  "TGFBI",
                  "VIM-PW", "NCAM1-PW", "SNAI1", "PRRX1") %>% rev()
  label <- "post-epipeelout2"
}

{ # Top markers for the Epi peelout
markers <- read_xlsx("results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/custom_res1_markers.xlsx")
markers$cluster <- factor(markers$cluster, levels = names(cols))
markers$score <- as.numeric(markers$score)
markers$score4 <- as.numeric(markers$score4)

markers$cluster %>% unique()
View(markers)

{ # Extended Data Fig. 8a (Heatmap of 50 markers)
markerdf2plot <- markers %>%
  arrange(cluster, -(abs(avg_log2FC))) %>%
  group_by(cluster) %>%
  slice_head(n = 50) %>%
  group_by(gene) %>%
  summarise(
    Cluster = if_else(n_distinct(cluster) > 1, "Mix", first(cluster)),
    .groups = "drop") %>%
    arrange(factor(Cluster, levels = c(names(cols), "Mix")))


View(markerdf2plot)
library(tidyr)
genes2plot <- markerdf2plot$gene
gene_annot <- markerdf2plot %>% column_to_rownames("gene")

label <- "50markers_bylog2FC_SCALE.DATA"
}

{# Extended Data Fig. 6c (Dotplot, all epi clusters)
genes2plot <- markers %>%
  group_by(cluster) %>%
  arrange(cluster, -score, , -(abs(avg_log2FC))) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

label <- "5markers_byscore_thenlog2FC"
}

{ # Markers for mesenchymal EPDCs (Fig. 6a)
genes2plot <- markers %>%
  dplyr::filter(cluster %in% c("EPDC_A", "09", "10", "02", "11", "16")) %>%
  group_by(cluster) %>%
  arrange(cluster, -score, , -(abs(avg_log2FC))) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

# Add RNAscope markers
genes2plot <- append(genes2plot, "FBLN1", after = 15) #FBLN1 after TGM2
genes2plot <- append(genes2plot, "PHEX", after = 3) #PHEX after LIX1
genes2plot
}

{ # Extended Data Fig. 4c, Extended Data Fig. 5a: QC plots
genes2plot <- c("S.Score", "G2M.Score",
                "nFeature_RNA", "nCount_RNA",
                "scrublet_Scores")
}

{ # Extended Data Fig.4d Big UMAP canonical marker dotplot
genes2plot <- c("TTN", "TNNT2", "MYH7-LOC138299692", "HEY2",
               "MYL7", "MYH6-LOC138299694",
               "NPPA-LOC138301123", "NPPA-PW", "HCN4",
               "HBA1", "HBA2",
               "CSF1R", "AIF1", "C1QA", "C1QB", "C1QC", "CD68",
               "MRC1-LOC138292840", "CD163-LOC138287444", "ANGPT2-LOC138284726", "TRPC3", "PLEK", 
               "CD79A", "CD79B", "PAX5", "CD40",
               "SKAP1", "CAMK4", "IL7R",
               "WT1", "TBX18", "TCF21", "BNC1",
               "LRRN4", "PODXL", "CLDN15", "CLDN6-LOC138304364", "CLDN7",
               "VIM-PW", "NCAM1-PW", "SNAI1", "PRRX1",
               "RGS5-LOC138293112", "NRXN1", "KCNJ8", "NRXN2",
               "MFAP4", "PRDM6", "ELN-PW",
               "MYH11", "ACTA2", "TAGLN",
               "VWA2", "SOX9",
               "POSTN", "DPT", "DCN",
               "COL1A1-PW", "COL1A2-PW", "COL1A2-LOC138268367",
               "COL3A1-PW", "COL5A1", "COL5A2", "COL5A3",
               "COL12A1", "VIM-PW", "FN1", "DKK2",
               "APLN", "KDR", "SOX18", "FLT1", "PECAM1",
               "CDH5", "VWF-PW", "LYVE1", "NFATC1", "NPR3", "TEK",
               "SOX10", "MPZ-PW", "MAG-LOC138246415", "NGFR-PW") %>% rev()
label <- "selectedmarkers"
}

{ # claudins # Extended Data Fig. 6
# rownames(obj2)[grep("CLDN", rownames(obj2))] %>% sort() %>% View()
genes2plot <- c("CLDN1", "CLDN3-LOC138286503",
"CLDN3-LOC138304359", "CLDN4-LOC138286500", "CLDN4-LOC138286501",
"CLDN4-LOC138304180", "CLDN5", "CLDN6-LOC138304362", "CLDN6-LOC138304364",
"CLDN7", "CLDN8-LOC138249872", "CLDN9-LOC138304358", "CLDN9-LOC138304361",
"CLDN10-LOC138249171", "CLDN10-LOC138265823", "CLDN10-LOC138303642",
"CLDN11", "CLDN12", "CLDN15", "CLDN16-LOC138302076", "CLDN16-LOC138302104",
"CLDN16-PW", "CLDN18", "CLDN19", "CLDN20", "CLDN23", "CLDN34")

label <- "claudins"
}

########################### Epicardial gene Module
{ # Extended Data Figure. 5e Epicardial lineage module to justify peelout
genes2plot <- c("WT1", "TBX18", "TCF21", "BNC1",
                "CLDN15", "CLDN6-LOC138304364", "CLDN7",
                "TGFBI", "PRRX1")
label <- "epi_lineage"
obj2 <- AddModuleScore(obj2, features = list(genes2plot), name = label)
genes2plot <- c("epi_lineage1", genes2plot) %>% rev()
}

########################## Basic quick feature plot
FeaturePlot(obj2,
            pt.size = ptsize,
            reduction = reduction,
            features = "epicardial_lineagev1",
            split.by = "orig.ident",
            order = TRUE) &
    scale_colour_viridis(option = "viridis")
# 
ggsave("results/figures/epi_lineage.png", width = 28)

########################### Feature Plots for list

# One FeaturePlot per gene
for (g in genes2plot) {
  p1 <- FeaturePlot(obj2,
                    pt.size = ptsize,
                    reduction = reduction,
                    features = g,
                    order = TRUE) +
    scale_colour_viridis(option = "viridis")

  # PNG save
  ggsave(filename = paste0("results/figures/", datasetname, "_", g, ".png"),
         width = 10, height = 10,
         plot = p1)
}

################################ Big UMAP Dotplot
# obj3 <- subset(obj2, named_clusters %in% c("EPDC_A", "09", "10", "02", "11", "16"))
p <- DotPlot(obj2,
        scale = TRUE,
        unique(genes2plot) %>% rev(),
        # idents = names(cols) %>% rev(),
        group.by = clustering) +
  scale_colour_viridis(option = "viridis") +
  scale_x_discrete(labels = function(x) stringr::str_to_sentence(x)) +
  coord_flip() +
  RotatedAxis()
p

ggsave(paste0("results/figures/", label, "_", datasetname, "_dotplotlot.pdf"),
       height = 13,
       width = 7,
       units = "in")

ggsave(paste0("results/figures/", label, "_", datasetname, "flipdotplot.pdf"),
       height = 0.15 * length(genes2plot) + 3,
       width = 7,
       units = "in")

########################### HEATMAP prep: list of genes
# Define color palette
my_palette <- colorRampPalette(plasma(20))(100)

# Scale all genes in the Seurat object (the original only scaled Variable Features)
obj3 <- DietSeurat(obj2, layers = c("counts", "data"))
obj3 <- ScaleData(obj3, features = rownames(obj3))

genes2plot <- rownames(gene_annot)[gene_annot$Cluster %in% c("EPDC_A", "09", "10", "02")]

# Specify order of cells: by clusters and then by condition
cellsordered <- obj2@meta.data %>%
  # dplyr::filter(!!sym(clustering) %in% c("03", "07", "04", "13", "14", "15", "12")) %>%
  # dplyr::filter(!(!!sym(clustering) %in% c("15", "12"))) %>%
  dplyr::filter(!!sym(clustering) %in% c("EPDC_A", "02", "09", "10")) %>%
  # dplyr::filter(!(orig.ident %in% c("10dpci_cCpE"))) %>%
  arrange(!!sym(clustering), orig.ident) %>%
  arrange(orig.ident, !!sym(clustering)) %>%
  rownames()

# Get data for custom heatmap
mat <- FetchData(obj3, genes2plot, slot = "scale.data") %>% t() # Expression data
mat_annot <- FetchData(obj3, c("orig.ident", clustering)) # Annotation data
colnames(mat_annot) <- c("Condition", "Cluster") # Identifying annotation data
View(mat_annot)
View(mat)

# Define colors for annotations
names(cols)
ann_colors <- list(Condition = c("PM_naive" = "#b7fff9",
                                 "14dpsham" = "#ffacf8",
                                 "10dpci_veh" = "#b2ffab",
                                 "10dpci_cCpE" = "#fffbc0",
                                 "14dpci" = "#cbabff"),
                   Cluster = c(cols, "Mix" = "brown"),
                   avg_log2FC = colorRampPalette(c("blue", "white", "red"))(100),
                   p_val_adj  = colorRampPalette(c(brewer.bupu(10)))(100))

# Choose: gap by clustering or by condition
# mat_annot$Cluster_condition <- paste0(mat_annot$Cluster, mat_annot$Condition)
mygaps <- which(!duplicated(mat_annot[cellsordered, ]$Cluster, fromLast = TRUE))
mygaps <- which(!duplicated(mat_annot[cellsordered, ]$Condition, fromLast = TRUE))
mygaps <- which(!duplicated(mat_annot[cellsordered, ]$Cluster_condition, fromLast = TRUE))

filename <- paste0("results/figures/", datasetname, "CHANGETHIS.pdf")
rownames(mat) <- stringr::str_to_sentence(rownames(mat))
rownames(gene_annot) <- stringr::str_to_sentence(rownames(gene_annot))

dim(mat)

# Make the heatmap
pheatmap(mat[, cellsordered],
         color = my_palette,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         show_rownames = FALSE,
         show_colnames = FALSE,
         angle_col = 45,
         gaps_col = mygaps,
        #  gaps_row = rowgaps,
         treeheight_row = 0,
         annotation_row = gene_annot,
         annotation_col = mat_annot,
         annotation_colors = ann_colors,
         cellheight = 1,
         cellwidth = 0.1,
         main = label,
        #  width = 8,
        #  height = 7,
         filename = filename
         )

dev.off()

############################### Violin Plots
# For any umap
VlnPlot(obj2, genes2plot, same.y.lims = TRUE, cols = cols, alpha = 0, ncol = 5) +
  RotatedAxis() &
  labs(x = "", y = "")
ggsave(paste0("results/figures/", datasetname, label, "nopts_vlnplot.png"), height = 12, width = 30)

######################### Proportions bar chart
# for epi only

# Fetch & format info on proportions
m <- table(obj2@meta.data[[clustering]], obj2$orig.ident)
m <- m[, colSums(m) > 0] # removing empty info about devo & 14dp conditions
df <- as.matrix(m) %>% as.data.frame(.)
colnames(df) <- c("Cluster", "Condition", "Frequency")

df$Cluster <- factor(df$Cluster, levels = names(cols))

df %>% View()

###########
# Make barplot (all clusters). Fig, 5f and Extended Data Fig. 5b
ggplot(df, aes(x = Condition, y = Frequency, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::percent_format())  +
  labs(y = "Percentage",
       title = datasetname) +
  # coord_flip() +
  theme_classic() +
  RotatedAxis()

ggsave(paste0("results/figures/", datasetname, "cluster_composition.png"),
       height = 10, width = 6)

###########
# Make piechart (mesothelial epi only). Fig. 5g
epiclusters <- c("03", "07", "04", "13", "14")
dfepi <- df %>% filter(Cluster %in% epiclusters)

for (condition in unique(dfepi$Condition)) {
  percon <- dfepi %>% dplyr::filter(Condition == condition)
  
  pdf(paste0("results/figures/10dpci_Epipeelout_", condition, "_Epipie.pdf"))
  pie(percon$Frequency, labels = percon$Cluster, col = cols[percon$Cluster])
  dev.off()
}

###########
# Make pie charts (EPDCs only). Fig. 6b
EPDCclusters <- c("02", "10", "09", "11", "16", "EPDC_A")

dfEPDC <- df %>% filter(Cluster %in% EPDCclusters)

for (condition in unique(dfEPDC$Condition)) {
  percon <- dfEPDC %>% dplyr::filter(Condition == condition)
  
  pdf(paste0("results/figures/10dpci_Epipeelout_", condition, "_EPDCpie.pdf"))
  pie(percon$Frequency, labels = percon$Cluster, col = cols[percon$Cluster])
  dev.off()
}

################################ Fig. 6d EPDC-A sham vs. cryoinjury (veh) ###########
my_palette <- colorRampPalette(viridis(20))(100) # custom viridis colors
options(enrichplot.colours = my_palette) # set viridis as colorscale

# Read in enrichment object
gomaster <- readRDS("results/DGE/ALGORAB/ALGORABEPDC_A_DGE_results_enrichment.rds")

# File with terms of interest
filename <- "results/DGE/ALGORAB/ALGORABEPDC_A_DGE_results_enrichment_annotatedbyInika.xlsx"

{ # Upregulated terms
go_obj_up <- gomaster[[3]]
sheet_up <- "GO_ORA_up"
label_up <- "ALGORAB_EPDCA_up"

# Get terms of interest
terms2plot_up <- read_xlsx(filename, sheet = sheet_up) %>%
  dplyr::filter(Keep == "Yes") %>%
  dplyr::pull(Description)

# Plot terms of interest
p_up_og <- dotplot(go_obj_up, showCategory = terms2plot_up)
# ggsave(paste0("results/figures/", label, ".pdf"), height = 5, width = 7)

}

{ # Downregulated terms
go_obj_down <- gomaster[[4]]
sheet_down <- "GO_ORA_down"
label_down <- "ALGORAB_EPDCA_down"

# Get terms of interest
terms2plot_down <- read_xlsx(filename, sheet = sheet_down) %>%
  dplyr::filter(Keep == "Yes") %>%
  dplyr::pull(Description)

# Plot terms of interest
p_down_og <- dotplot(go_obj_down, showCategory = terms2plot_down)
# ggsave(paste0("results/figures/", label, ".pdf"), height = 5, width = 7)

}

p_up_og + p_down_og

# Extract results
df_up <- gomaster[[3]]@result %>%
  filter(Description %in% terms2plot_up) %>%
  mutate(Direction = "Upregulated")

df_down <- gomaster[[4]]@result %>%
  filter(Description %in% terms2plot_down) %>%
  mutate(Direction = "Downregulated")

# Combine
df <- bind_rows(df_up, df_down) %>%
  # Parse GeneRatio (e.g. "10/203")
  tidyr::separate(GeneRatio, into = c("hits", "total"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio = hits / total) %>%
  select(-hits, -total) %>%
  group_by(Direction) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  # mutate(Description = reorder_within(Description, GeneRatio, Direction)) %>%
  ungroup()


order <- df %>% arrange(GeneRatio) %>% dplyr::pull(Description) %>% unique()
df$Description <- factor(df$Description, levels = order)
df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))

# Plot
ggplot(df,
       aes(x = GeneRatio,
           y = Description,
           colour = p.adjust,
           size = Count)) +
  geom_point() +
  facet_wrap(~Direction, ncol = 1, scales = "free_y") +
  scale_colour_gradientn(
    colours = my_palette,
    trans = "reverse",      # optional: small p-values at top of legend
    breaks = c(0.001, 0.003, 0.005, 0.007)
  ) +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = NULL,
    colour = "Adjusted p-value"
  ) +
  scale_y_discrete(position = "right")

ggsave("results/figures/epdc_a_dotplot.pdf", height = 3.5, width = 8)

################# Fig. 6e: Sham vs Veh: Collect genes for term of interest
toi <- "extracellular structure organization"
goi <- df[df$Description == toi,]$geneID

genes <- strsplit(goi, "/") %>% unlist()

# Load unranked genelist to fetch the DEGs
identifiers_spec <- read_xlsx(filename, sheet = "Unranked_sig_genelist") %>%
  dplyr::select(c("genename_Hs", "gene"))

genes_Pw <- identifiers_spec %>%
  dplyr::filter(genename_Hs %in% genes) %>%
  dplyr::pull(gene)

genes2plot <- genes_Pw

# Specify colors for the violin plot
cond_cols <- c("PM_naive" = "#b7fff9",
               "14dpsham" = "#ffacf8",
               "10dpci_veh" = "#b2ffab",
               "10dpci_cCpE" = "#fffbc0",
               "14dpci" = "#cbabff")

# Subset object to only include sham and vehicle
obj3 <- subset(obj2, subset = orig.ident %in% c("14dpsham", "10dpci_veh"))

# Break up the list of top_genes into 20-gene chunks so printing is easy
gene_groups <- split(genes2plot, ceiling(seq_along(genes2plot) / 20))
pdf("results/figures/ECM_EPDC_A_violinplots_ALGORABgenes.pdf", height = 10, width = 10)
for (i in gene_groups) {
p <- VlnPlot(obj3,
             #  cols = cond_clust_cols,
             features = i,
             ncol = 5,
             alpha = 0,
             cols = cond_cols,
             group.by = "orig.ident",
             idents = c("EPDC_A")
             ) &
            #  theme_void() &
        theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())
  print(p)
}

dev.off()

##################### Pseudobulk expression (Fig. 7a RNAscope heatmap)
# Identify cluster & condition groups with <3 cells
freq_tab <- table(obj2$orig.ident, obj2$merged_clusters) %>%
  as.data.frame() %>%
  rename(., Condition = Var1) %>%
  rename(., Cluster = Var2) %>%
  mutate(Condition_Cluster = paste0(Condition, "_", Cluster))

torem <- freq_tab %>% dplyr::filter(Freq < 3)
torem

# Remove groups with < 3 cells
obj2$Condition_Cluster <- paste0(obj2$orig.ident, "_", obj2$merged_clusters)
obj3 <- subset(obj2,
               subset = !(Condition_Cluster %in% torem$Condition_Cluster))

DimPlot(obj3, group.by = "orig.ident")


# Aggregare the expression
pb <- AggregateExpression(obj3,
                          assays = "RNA",
                          group.by = c("orig.ident", "merged_clusters"),
                          return.seurat = TRUE)

pb@meta.data %>% View()

# Re-add the orig.ident column
pb@meta.data$pb_idents <- pb@meta.data$orig.ident # Save the pb idents
pb@meta.data$orig.ident <- sub("_.*", "", pb@meta.data$orig.ident) #Remove allt efter first underscore
pb@meta.data$orig.ident <- gsub("g", "", pb@meta.data$orig.ident) #Remove the "g" preceding 10dpci-veh and 14dpci
unique(pb@meta.data$orig.ident)
pb@meta.data$orig.ident <- factor(pb@meta.data$orig.ident,
                                  levels = c("PM-naive",
                                             "14dpsham",
                                             "10dpci-veh",
                                             "10dpci-cCpE"))

# Make named clusters match the known levels
pb@meta.data$merged_clusters <- gsub("-", "_", pb@meta.data$merged_clusters)
pb@meta.data$merged_clusters <- factor(pb@meta.data$merged_clusters,
                                      levels = unname(map))

# Define colors for annotations
names(cols)
ann_colors <- list(Condition = c("PM-naive" = "#b7fff9",
                                 "14dpsham" = "#ffacf8",
                                 "10dpci-veh" = "#b2ffab",
                                 "10dpci-cCpE" = "#fffbc0",
                                 "14dpci" = "#cbabff"),
                   Cluster = cols)

# Specify palette
my_palette <- colorRampPalette(plasma(20))(100)

# Specify order of cells:
cellsordered <- pb@meta.data %>%
  # dplyr::filter(!(!!sym(clustering) %in% c("03", "07", "04"))) %>%
  arrange(orig.ident, merged_clusters) %>% #by condition & then by cluster
  arrange(merged_clusters, orig.ident) %>% #by condition & then by cluster
  rownames()

# Select genes to plot
genes2plot <- top_genes
genes2plot <- c("BNC1", "PHEX", "CILP", "TNC-PW", "FBLN1") %>% rev()

{ # Data slot
mat <- FetchData(pb, genes2plot %>% rev(), slot = "data") %>% t() # Log-normalized Expression Data
mytitle <- "Data slot"
}

# Specify annotations
mat_annot <- FetchData(pb, c("orig.ident", "merged_clusters")) # Annotation data
colnames(mat_annot) <- c("Condition", "Cluster")
head(mat_annot)

# Choose: gap by clustering or by condition
mygaps <- which(!duplicated(mat_annot[cellsordered, ]$Cluster, fromLast = TRUE))
mygaps <- which(!duplicated(mat_annot[cellsordered, ]$Condition, fromLast = TRUE))

# Make the gene names match convention
rownames(mat) <- stringr::str_to_sentence(rownames(mat))

# Specify filename to save heatmap to
filename <- paste0("results/figures/", datasetname, "CHANGETHIS4.pdf")

# Make the heatmap
pheatmap(mat[, cellsordered],
         color = my_palette,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
        #  show_rownames = FALSE,
         show_colnames = TRUE,
         main = mytitle,
         angle_col = 45,
         gaps_col = mygaps,
         treeheight_row = 0,
         scale = "none",
         annotation_col = mat_annot,
         annotation_colors = ann_colors,
         cellheight = 12,
         cellwidth = 12,
        #  width = 8,
        #  height = 5,
         filename = filename)

dev.off()
