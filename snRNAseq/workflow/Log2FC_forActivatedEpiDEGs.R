##################### README
# Author: Inika Prasad, inika.prasad@gmail.com
# Date: 2026-04-26
# Aim: we performed Differential Gene Expression on activated epicardium (cluster 04, which is called cluster 03 in the paper)
# We wanted to visualize log2FC values significant veh-ccpe DEGS in both sham-veh and veh-ccpe.
# But because some genes were expressed in too few genes to be considered for the DGE analysis,
# They have NA values in the log2FC column of the sham-veh comparison.

###################### Load libraries ################################
library(tibble)

###################### Prep work
# Run DEG_heatmaps_prep.R for all_degs dataframe
# (to check whether the indpendent calculation of log2FC below matches)

####################### Load data ####################################

vehccpe_degs_cluster04 <- read_xlsx("results/DGE/ALGOL/ALGOL04_DGE_results.xlsx", sheet = "sig_degs")
vehccpe_degs <- vehccpe_degs_cluster04$gene

input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
clustering <- "merged_clusters"

obj <- readRDS(input_rds)
obj2 <- JoinLayers(obj)
Idents(obj2) <- clustering
obj2@meta.data[[clustering]] <- factor(obj2@meta.data[[clustering]], levels = names(cols))
levels(obj2)

DimPlot(obj2)

# Make metadata for condition & cluster
obj2$condition_cluster <- paste0(obj2$orig.ident, obj2$merged_clusters)
Idents(obj2) <- "condition_cluster"

###################### Compute log2 Fold Changes anew #################


# Find veh vs. ccpe (ALGOL) fold changes
ALGOL_attempt <- FoldChange(object = obj2, ident.1 = "10dpci_cCpE04", ident.2 = "10dpci_veh04", features = vehccpe_degs)
colnames(ALGOL_attempt) <- paste0(colnames(ALGOL_attempt), "_ALGOLattempt")
ALGOL_attempt <- ALGOL_attempt %>% rownames_to_column("gene")

df <- left_join(all_degs, ALGOL_attempt, by = "gene")

# Find sham vs. vehicle (ALGORAB) fold changes
ALGORAB_attempt <- FoldChange(object = obj2, ident.1 = "10dpci_veh04", ident.2 = "14dpsham04", features = vehccpe_degs)
colnames(ALGORAB_attempt) <- paste0(colnames(ALGORAB_attempt), "_ALGORABattempt")
ALGORAB_attempt <- ALGORAB_attempt %>% rownames_to_column("gene")
View(ALGORAB_attempt)

df2 <- left_join(df, ALGORAB_attempt, by = "gene")

df2 %>% View()
df2$ALGORAB04 == df2$avg_log2FC_ALGORABattempt
df2$ALGORAB04[37] == df2$avg_log2FC_ALGORABattempt[37]

##################### Rewrite all_degs object for use in heatmap ###############
all_degs <- df2 %>%
    dplyr::select(c(gene, avg_log2FC_ALGORABattempt, avg_log2FC_ALGOLattempt))

all_degs %>% colnames()
colnames(all_degs) <- c("gene", "ALGORAB", "ALGOL")

###############
# To make the heatmap, go to DEG_heatmaps.R