##################### README
# Author: Inika Prasad, inika.prasad@gmail.com
# Date: 2026-04-26
# Aim: we performed Differential Gene Expression on activated epicardium (cluster 04, which is called cluster 03 in the paper)
# We wanted to visualize log2FC values significant veh-ccpe DEGS in both sham-veh and veh-ccpe.
# But because some genes were expressed in too few genes to be considered for the DGE analysis,
# They have NA values in the log2FC column of the sham-veh comparison.

###################### Load libraries ################################
library(tibble)
library(stringr)
library(dplyr)
library(org.Hs.eg.db)
library(readxl)
library(Seurat)
library(pheatmap)
library(plotly)
library(scales)
library(writexl)

set.seed(100)

###################### Prep work
# Run DEG_heatmaps_prep.R for all_degs dataframe
# (to check whether the indpendent calculation of log2FC below matches)

####################### Load data ####################################
# Load existing info about veh vs. cCpE DEGs
vehccpe_degs_cluster04 <- read_xlsx("results/DGE/ALGOL/ALGOL04_DGE_results.xlsx", sheet = "sig_degs")
vehccpe_degs <- vehccpe_degs_cluster04$gene

# Specify RDS object with Epi Peelout
input_rds <-"results/16_peeled_integrated_DimRed/10dpci_ctrl/10dpci_ctrl_DRscrub_Epi_dimred_CUSTOMRES.rds"
clustering <- "merged_clusters"

# Specify metadata

# Load RDS object
obj <- readRDS(input_rds)
obj2 <- JoinLayers(obj)
Idents(obj2) <- clustering
levels(obj2)

# Vibe check the RDS object
DimPlot(obj2, label = TRUE)

# Make metadata for condition & cluster
obj2$condition_cluster <- paste0(obj2$orig.ident, obj2$merged_clusters)
Idents(obj2) <- "condition_cluster"
levels(obj2)

###################### Compute log2 Fold Changes anew #################

# Find veh vs. ccpe (ALGOL) fold changes
ALGOL_attempt <- FoldChange(object = obj2,
                            ident.1 = "10dpci_cCpE04",
                            ident.2 = "10dpci_veh04",
                            features = vehccpe_degs)
colnames(ALGOL_attempt) <- paste0(colnames(ALGOL_attempt), "_ALGOLattempt")
ALGOL_attempt <- ALGOL_attempt %>% rownames_to_column("gene")

# Combine de novo computed log2FC with existing pre-calculated ones
df <- left_join(vehccpe_degs_cluster04, ALGOL_attempt, by = "gene")

# Find sham vs. vehicle (ALGORAB) fold changes
ALGORAB_attempt <- FoldChange(object = obj2,
                              ident.1 = "10dpci_veh04",
                              ident.2 = "14dpsham04",
                              features = vehccpe_degs)
colnames(ALGORAB_attempt) <- paste0(colnames(ALGORAB_attempt), "_ALGORABattempt")
ALGORAB_attempt <- ALGORAB_attempt %>% rownames_to_column("gene")
View(ALGORAB_attempt)

# Cobine veh vs. ccpe with sham vs. veh
df2 <- left_join(df, ALGORAB_attempt, by = "gene")

# Sanity check: old and new computed fold changes match for veh-cCpE
# They almost match, there is a miniscule error (10^-16, which disappears after rounding)
df2 %>% View()
df2 %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    dplyr::filter(avg_log2FC != avg_log2FC_ALGOLattempt) %>%
    mutate(checkdiff = avg_log2FC - avg_log2FC_ALGOLattempt) %>%
    View() # Should be empty

################################ Load identifiers
identifiers <- read_xlsx("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.xlsx")
identifiermap <- setNames(identifiers$gene_id, identifiers$genename_Hs_PwNameFormat)


# Round off the values so that the match is exact
# And convert gene names back to the original NCBI
df3 <- df2 %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(gene_id = identifiermap[gene])

################################ Make a heatmap

mat_data <- df3 %>%
    dplyr::select(avg_log2FC_ALGOLattempt, avg_log2FC_ALGORABattempt) %>%
    as.matrix()

# Make gene names Sentence case, and the right identifiers
rownames(mat_data) <- str_to_sentence(df3$gene_id)

# Specify filename
myfilename <- "results/figures/ActivatedEpi_vehccpeDEGs.pdf"

# Specify scale end
scale_end <- max(abs(max(mat_data)), abs(min(mat_data))) %>% ceiling()

# Make heatmap
pheatmap(t(mat_data),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 10,
         fontsize_col = 14,
         angle_col = 45,
         gaps_col = mygaps,
         show_rownames = TRUE,
         treeheight_row = 0,
         treeheight_col = 0,
         border_color = "white",
         fontsize = 7,
         cellheight = 30,
         cellwidth = 22,
         legend_breaks = c(-scale_end, -scale_end/2, 0, scale_end/2, scale_end),
         width = 23,
         height = 8,
         filename = myfilename,
        )

################################ Make a sankey plot
# Find GO Term IDs associated with genes
go_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = df3$gene_id,
                                columns = c("GO"),
                                keytype = "SYMBOL")
unique(go_map$SYMBOL)

# Remove NAs arising from unmapped genes
go_map <- go_map[!is.na(go_map$GO), ]
go_map %>% View()

# Find human-readable descriptions for the GO IDs
go_map$GO_Description <- Term(go_map$GO)
View(go_map)

# Which terms have many genes associated?
table(go_map$GO_Description)[table(go_map$GO_Description) > 1] %>% sort()
unique(go_map$SYMBOL) # 41 genes

############################### Make Supplementary file
# Rename columns of df3
new_names <- setNames(colnames(df3)[!grepl("score|ALGOLattempt", colnames(df3))],
  colnames(df3) %>%
    .[!grepl("score|ALGOLattempt", .)] %>% #remove score & redundant columns
    gsub("ALGORABattempt", "Sham_Veh", .) %>% # Replace ALGORAB with Sham_Veh
    ifelse(
      . %in% c("gene", "gene_id") | grepl("Sham_Veh", .),
      .,
      paste0(., "_Veh_cCpE") # Add Veh_cCpE to ALGOL info
    ))

df4 <- df3 %>%
  dplyr::select(-contains("score"), -contains("ALGOLattempt")) %>%
  rename(!!!new_names) %>%
  relocate(gene_id, gene)

View(df3)
View(df4)

write_xlsx(list(DEGs_ActivatedEpi = df4,
                AssociatedGOterms = go_map),
           path = "results/figures/ActivatedEpi_DEGs_Supplementary.xlsx")


################################ Identify genes and processes of interest

# OPTION A: Specify genes of interest
goi <- c("COL24A1", "KRT8", "OMA1", "FABP7", "CD36", "PTGS2",
         "EPHA2", "CRISPLD2", "STYX", "RASGEF1A", "SYNPO", "PLEKHG2",
         "EPHA10", "ERG", "PTPN14", "HR", "CRIP2", "WT1")

go_map_subset <- go_map %>%
    dplyr::filter(SYMBOL %in% goi)

# OPTION B: Don't specify genes of interest
go_map_subset <- go_map


# Define how GO terms should be grouped
keywords <- list(
  Ephrin_signaling = c("ephrin"),
  Erk_12_signaling_pathway = c("ERK"),
  Transcriptional_regulation = c("regulation of transcription"),
  Cyto_migration = c("cytoskeleton", "migration"),
  Junction_adh_ecm = c("junction", "adhesion", "extracellular matrix"),
  Lipid_metabolism = c("lipid"),
  Cell_differentiation = c("differentiation"))

# Map from GO term to GO category
go_df <- go_map_subset %>%
  rowwise() %>%
  mutate(Category = names(keywords)[
      which(sapply(keywords, \(k) str_detect(GO_Description, paste(k, collapse = "|"))))[1]
    ]) %>%
  ungroup()

View(go_df)

# Add custom categories
custom_cat <- c(PTPN14 = "Hippo_yap", # Liu et al. 2012, Nature Oncogene
                CRISPLD2 = "Hippo_yap", # Zhang et al, 2025, Int J Biol Sci
                ERG = "Erk_12_signaling_pathway", # Kedage et al, 2017, J Biol Chem
                CRIP2 = "Transcriptional_regulation", #Cheung et al, 2011, 10.1073/pnas.1101747108
                PTGS2 = "Lipid_metabolism",
                RASGEF1A = "Erk_12_signaling_pathway" # NEEDED
)

lit_df <- data.frame(SYMBOL = names(custom_cat),
           Category = unname(custom_cat))

# Combine GO term & literature categories
df <- bind_rows(go_df, lit_df) %>%
  mutate(SYMBOL_og = SYMBOL,
         SYMBOL = str_to_sentence(SYMBOL))

####################### Make a Sankey plot

edges <- df %>%
    dplyr::select(sourceDescription = SYMBOL,
                  targetDescription = Category,
                  SYMBOL_og) %>%
    distinct() %>%
    dplyr::filter(!is.na(targetDescription))

# View(edges)

# READABLE, Specify plot details with GO term names
nodes <- data.frame(name = unique(c(edges$SYMBOL_og, edges$targetDescription))) # USE original SYMBOLS so that merge works
nodes <- left_join(nodes, df3, by = join_by(name == gene_id))
nodes <- nodes %>%
  mutate(name_og = name,
         name = str_to_sentence(name))
# nodes$avg_log2FC[is.na(nodes$avg_log2FC)] <- 0
View(nodes)

# Order Left-hand nodes
left_nodes <- c("Fabp7", "Oma1", "St6galnac5", "Rnf213",
  "Abca1", "Hr", "Hif3a", "Crip2", "Erg", "Wt1", "Slc25a5",
  "Nbr1", "Ptgs2", "Cd36", "Gapdh", "F3", "Krt8",
  "Ppfia2", "Col24a1", "Atp2a2", "Epha2", "Synpo",
  "Rasgef1a", "Crispld2", "Styx", "Pde8a", "Ptpn14",
  "Epha10")

# Order Right-hand nodes
right_nodes <- c("Lipid_metabolism", "Transcriptional_regulation",
  "Cell_differentiation",  "Junction_adh_ecm", "Cyto_migration",
  "Erk_12_signaling_pathway", "Ephrin_signaling", "Hippo_yap")

myorder <- c(left_nodes, right_nodes)
setdiff(nodes$name, myorder) # should be 0
setdiff(myorder, nodes$name) # should be 0

nodes <- nodes %>%
  dplyr::arrange(match(name, myorder)) %>%
  mutate(x = c(rep(0.05, length(left_nodes)), rep(0.90, length(right_nodes))),
         y = c(seq(0.02, 0.98, length.out = length(left_nodes)),
               seq(0.08, 0.99, length.out = length(right_nodes))))


# Create node IDs
edges$source_id <- match(edges$sourceDescription, nodes$name) - 1
edges$target_id <- match(edges$targetDescription, nodes$name) - 1
lim <- max(abs(nodes$avg_log2FC))

# Make node size uniform or variable
edges <- edges %>%
  group_by(sourceDescription) %>%
#   mutate(nlinks = 1 / n()) %>% # to make nodes uniform
  mutate(nlinks = 1) %>% # to make nodes variable
  ungroup()

# Specify edge colour based on target node
pal2 <- col_factor(
  palette = pal_hue(l = 50)(length(unique(edges$targetDescription))),
  domain = as.factor(edges$targetDescription),
  reverse = TRUE
)

edges$color <- pal2(as.factor(edges$targetDescription))
edges$color <- alpha(edges$color, 0.7) # Specify transparency for the colors
show_col(edges$color)

# Color the Gene nodes by log2FC value
pal <- col_numeric(
  palette = "RdYlBu",
  domain = c(min(mat_data), max(mat_data)), # So it matches the heatmap
  reverse = TRUE,
  na.color = NA) # The pathway nodes will be white)

nodes$color <- pal(nodes$avg_log2FC)
show_col(nodes$color)
View(nodes)

# Color the pathway nodes by incoming link color where possible
pathway_colmap <- edges %>%
  dplyr::select(targetDescription, color) %>%
  distinct()

nodes <- nodes %>%
  left_join(pathway_colmap, by = join_by(name == targetDescription), suffix = c("", "_new")) %>%
  mutate(color = coalesce(color_new, color)) %>%
  select(-color_new)

View(nodes)

# Build plot
sankeyp <- plot_ly(
  type = "sankey",
  node = list(
    label = nodes$name,
    # label = rep("", length(nodes$name)),
    color = nodes$color,
    height = 8,
    x = nodes$x,
    y = nodes$y,
    customdata = nodes$avg_log2FC,
    hovertemplate = paste0(
      "%{label}<br>",
      "log2FC: %{customdata:.3f}<extra></extra>")),
  link = list(
    color = edges$color,
    source = edges$source_id,
    target = edges$target_id,
    value = edges$nlinks)
)

sankeyp 

library(reticulate)
use_condaenv("snRNAseq_2025_env", required = TRUE)
py_config()
plotly::save_image(sankeyp, "results/figures/cluster04_DEGs_sankey2.pdf",
  width = 300, height = 1000)
