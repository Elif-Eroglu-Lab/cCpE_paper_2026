################################ README
# Figures for veh-ccpe manuscript
# Author: Inika Prasad
# Last modified: 2026-09-28

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

map <- setNames(df$Naming_manuscript, df$Cluster) # (desired naming, original naming)
map

cols <- setNames(df$cols, df$Naming_manuscript) # (colors, desired naming)

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


map <- setNames(df$Naming_manuscript, df$Cluster) # (desired naming, original naming)
map

cols <- setNames(df$cols, df$Naming_manuscript) # (colors, desired naming)

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

################################ Load identifiers
identifiers <- read_xlsx("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.xlsx")
identifiermap <- setNames(identifiers$gene_id, identifiers$genename_Hs_PwNameFormat)

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
  label <- "post-epipeelout"
  height <- 5
  width <- 6
}

{ # Top markers for the Epi peelout
markers <- read_xlsx("results/16b_clustermarkers_4peeledintdata/10dpci_ctrl/custom_res1/custom_res1_markers.xlsx")
markers$cluster <- factor(markers$cluster, levels = df$Cluster)
markers$score <- as.numeric(markers$score)

# Make supplementary file for marker genes
  
marker_Supplementary <- markers %>%
  dplyr::select(c(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, abundance_diff, score)) %>%
  mutate(gene_id = identifiermap[gene],
         Cluster_manuscriptname = map[cluster]) %>%
         dplyr::arrange(Cluster_manuscriptname) %>%
         relocate(gene_id, Cluster_manuscriptname, gene, cluster)

write_xlsx(marker_Supplementary, path = "results/figures/EpiMarkers_supplementary.xlsx")

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
    arrange(factor(Cluster, levels = c(df$Cluster, "Mix")))


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
  unique() %>%
  rev()

label <- "5markers_byscore_thenlog2FC"
width <- 6.5
height <- 12
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
genes2plot <- genes2plot %>% rev()

label <- "EPDCs_topmarkers"
height <- 7
width <- 5
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

{ # Extended Data Fig. 4c, Extended Data Fig. 5a: QC plots
genes2plot <- c("S.Score", "G2M.Score",
                "nFeature_RNA", "nCount_RNA",
                "scrublet_Scores")
}

{ # Extended Data Fig.4d Big UMAP canonical marker dotplot
genes2plot <- c("TTN", "TNNT2", 
               "MYH7-LOC138299692", # myosin-7
               "HEY2",
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
width <- 10
height <- 16
}


########################### Epicardial gene Module
{ # Extended Data Figure. 5e Epicardial lineage module to justify peelout
genes2plot <- c("WT1", "TBX18", "TCF21", "BNC1",
                "CLDN15", "CLDN6-LOC138304364", "CLDN7",
                "TGFBI", "PRRX1")
label <- "epi_lineage"
obj2 <- AddModuleScore(obj2, features = list(genes2plot), name = label)
genes2plot <- c("epi_lineage1", genes2plot) %>% rev()
width <- 10
height <- 3
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
obj3 <- subset(obj2, named_clusters %in% c("06", "07", "08", "09", "10", "11")) # For Fig. 6a

p <- DotPlot(obj2,
        scale = TRUE,
        unique(genes2plot),
        group.by = clustering) +
  scale_colour_viridis(option = "viridis") +
  scale_x_discrete(labels = function(x) stringr::str_to_sentence(identifiermap[x])) +
  coord_flip() +
  RotatedAxis()
p

ggsave(paste0("results/figures/", label, "_", datasetname, "_dotplot.pdf"),
       height = height,
       width = width,
       units = "in")

########################### HEATMAP prep: list of genes
# Define color palette
my_palette <- colorRampPalette(plasma(20))(100)

# Scale all genes in the Seurat object (the original only scaled Variable Features)
obj3 <- DietSeurat(obj2, layers = c("counts", "data"))
obj3 <- ScaleData(obj3, features = rownames(obj3), vars.to.regress = c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score"))

# genes2plot <- rownames(gene_annot)[gene_annot$Cluster %in% c("EPDC_A", "09", "10", "02")]

# Specify order of cells: by clusters and then by condition
cellsordered <- obj2@meta.data %>%
  # dplyr::filter(!!sym(clustering) %in% c("03", "07", "04", "13", "14", "15", "12")) %>%
  # dplyr::filter(!(!!sym(clustering) %in% c("15", "12"))) %>%
  # dplyr::filter(!!sym(clustering) %in% c("EPDC_A", "02", "09", "10")) %>%
  # dplyr::filter(!(orig.ident %in% c("10dpci_cCpE"))) %>%
  arrange(!!sym(clustering), orig.ident) %>%
  # arrange(orig.ident, !!sym(clustering)) %>%
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
rownames(mat) <- stringr::str_to_sentence(identifiermap[rownames(mat)])
rownames(gene_annot) <- stringr::str_to_sentence(rownames(gene_annot))

dim(mat)
rownames(mat)

# Make the heatmap
pheatmap(mat[, cellsordered],
         color = my_palette,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 10,
         show_rownames = TRUE,
         show_colnames = FALSE,
         angle_col = 45,
         gaps_col = mygaps,
        #  gaps_row = rowgaps,
         treeheight_row = 0,
        #  annotation_row = gene_annot,
         annotation_col = mat_annot,
         annotation_colors = ann_colors,
         cellheight = 12,
         cellwidth = 0.2,
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
