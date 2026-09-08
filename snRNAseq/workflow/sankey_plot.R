############### README
# Author: Inika Prasad
# Date created: 2026-03-17
# Aim: to look at the GO terms for a list of genes
# in a Sankey plot
# Here, I want to look at differentially expressed genes 
# between Vehicle and c-CpE

############### Load libraries
library(readxl)
library(writexl)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tibble)
library(plotly)
library(scales)

################ Select dataset
{
log2fc_cutoff <- 0.6
coi <- "EPDC_A"
degfile <- "results/DGE/ALGOL/ALGOLEPDC_A_DGE_results.xlsx"
}

{
log2fc_cutoff <- 0.6
coi <- "04"
degfile <- "results/DGE/ALGOL/ALGOL04_DGE_results.xlsx"
go_file <- gsub(".xlsx", "_associated_GOterms_fulltable.xlsx", degfile)
termfile <- gsub(".xlsx", "_associated_GOterms.xlsx", degfile)
umbrellafile <- gsub(".xlsx", "_modified_plotedges.xlsx", degfile)
}

###############  Load data
# Read in significant DEGs
markerdf <- read_xlsx(degfile,
                      sheet = "sig_degs") %>%
  dplyr::filter(abs(avg_log2FC) > log2fc_cutoff) %>% # Add log2FC cutoff if needed
  arrange(desc(avg_log2FC))

######################### Map gene names to Human gene names
# Map identifiers to Human symbols
identifiers <- read.csv("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.tsv", sep = "\t")
colnames(identifiers)

identifiers <- identifiers %>% dplyr::select(c("genename_Hs_PwNameFormat", "genename_Hs"))

# Map human symbols to gene identifiers
markerdf2 <- left_join(markerdf, identifiers, join_by(gene == genename_Hs_PwNameFormat))
View(markerdf2)

################ Get GO pathways for genes of interest

go_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = markerdf2$genename_Hs,
                                columns = c("GO"),
                                keytype = "SYMBOL")

unique(go_map$SYMBOL)

# Remove NAs arising from unmapped genes
go_map <- go_map[!is.na(go_map$GO), ]
go_map %>% View()

###############  Get the hierarchy of pathways
go_parents <- as.list(GO.db::GOBPPARENTS)  # for Biological Process

get_go_parents <- function(go_ids) {
  edges <- data.frame()

  # Get the parent terms for the main term
  for (id in go_ids) {
    # print(id)
    parent_list <- go_parents[[id]]
    if (!is.null(parent_list)) {
    df <- data.frame(child = id,
                     parent = parent_list) %>% rownames_to_column("type")
    # print(df)
    edges <- bind_rows(edges, df)
    }

  }

  # Output the edges of the dataframe
  edges

}

##############  Build edges: gene to GO term edges
edges_gp <- go_map %>%
  left_join(markerdf, by = join_by(SYMBOL == gene)) %>%
  transmute(source = SYMBOL,
            target = GO)

edges_gp$targetDescription <- Term(edges_gp$target)
edges_gp$sourceDescription <- edges_gp$source # Placeholder, can be modified to Title case or original Pw gene names or something
head(edges_gp)

unique(edges_gp$source)

# Record GO Terms associated with each gene
write_xlsx(edges_gp, go_file)

############### Filter out and Manually add edges
# What pathways do we have? Are any of them particularly interesting
# Save the pathway list to select manually those that are of interest
table(edges_gp$targetDescription) %>%
  sort() %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("term")
  # %>%
  # write_xlsx(termfile)

# Get the pathways of interest and filter out everything else
terms2keep <- read_xlsx(termfile, sheet = "selectpathways") %>%
    filter(Retain == "Yes") %>%
    pull(term)
edges_gp <- edges_gp[edges_gp$targetDescription %in% terms2keep, ]

# Add some edges manually, see file for evidence of connections
manual_edges <- read_xlsx(termfile, sheet = "manual_connections") %>%
  dplyr::select(-c(notes, citation))
edges_gp <- rbind(edges_gp, manual_edges)

# Remove some nodes for cleaner message
nodes2remove <- read_xlsx(termfile, sheet = "remove_nodes") %>%
  pull(nodes2remove)
edges_gp <- edges_gp[!edges_gp$sourceDescription %in% nodes2remove, ]

############### Build edges: GO term to parent term
go_edges <- get_go_parents(unique(edges_gp$target))
View(go_edges)
edges_pp <- go_edges %>%
  transmute(source = child,
            target = parent)

edges_pp$targetDescription <- Term(edges_pp$target)
edges_pp$sourceDescription <- Term(edges_pp$source)
head(edges_pp)
table(edges_pp$targetDescription) %>% sort() %>% as.matrix() %>% View()

############### Modify nodes in the edges
# we want to group together things for aesthetic purposes
# This file was modified manually after saving
# DO NOT OVERWRITE

# write_xlsx(edges_gp %>%
#   mutate(targetOG = targetDescription) %>%
#   dplyr::select(c(targetOG, targetDescription)) %>%
#   distinct(),
#   umbrellafile)

################ Combine edges
# OPTION 1: Prep plot with genes and parent terms (Depth = 1)
edges <- edges_gp %>%
  distinct()

# OPTION 2: Prep plot with genes, parent terms, and grandparent terms (Depth = 2)
edges <- bind_rows(edges_gp, edges_pp) %>%
  distinct()

# OPTION 3: Prep plot with custom edges
term_groups <- read_xlsx(umbrellafile)
edges <- left_join(edges_gp,
                   term_groups,
                   suffix = c(".og", ""),
                   by = join_by(targetDescription == targetOG)) %>%
  # Remove target because it contains distinct GO terms
  # which makes otherwise identical rows distinct
  dplyr::select(-target) %>%
  distinct()

nrow(edges)
edges %>% View()

# ################# Filter edges
# # Only retain nodes that have 2 or more inputs?
# ELEPHANT this removes some children nodes
# ELEPHANT: NEEDS TO BE DONE AT THE PP and GP level
# tokeep <- table(edges$targetDescription)[table(edges$targetDescription) >= 2] %>%
#     as.matrix() %>%
#     as.data.frame() %>%
#     rownames()

# head(tokeep)

# edges <- edges[edges$targetDescription %in% tokeep, ]
# nrow(edges)

################# Build plot
# READABLE, Specify plot details with GO term names
nodes <- data.frame(name = unique(c(edges$sourceDescription, edges$targetDescription)))
nodes <- left_join(nodes, markerdf2, by = join_by(name == genename_Hs))
# nodes$avg_log2FC[is.na(nodes$avg_log2FC)] <- 0

# Create node IDs
edges$source_id <- match(edges$sourceDescription, nodes$name) - 1
edges$target_id <- match(edges$targetDescription, nodes$name) - 1
lim <- max(abs(nodes$avg_log2FC))

# Make node size uniform or variable
edges <- edges %>%
  group_by(sourceDescription) %>%
  # mutate(nlinks = 1 / n()) %>% # to make nodes uniform
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
  domain = nodes$avg_log2FC,
  reverse = TRUE,
  na.color = NA # The pathway nodes will be white
)
nodes$color <- pal(nodes$avg_log2FC)
show_col(nodes$color)
View(nodes)

# Color the pathway nodes by incoming link color where possible
pathway_colmap <- edges %>%
  dplyr::select(targetDescription, color) %>%
  distinct()
nodes %>% mutate(color2 = )

nodes <- nodes %>%
  left_join(pathway_colmap, by = join_by(name == targetDescription), suffix = c("", "_new")) %>%
  mutate(color = coalesce(color_new, color)) %>%
  select(-color_new)

# Build plot
sankeyp <- plot_ly(
  type = "sankey",
  node = list(
    label = nodes$name,
    # label = rep("", length(nodes$name)),
    color = nodes$color,
    height = 8,
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

library(reticulate)
use_condaenv("snRNAseq_2025_env", required = TRUE)
py_config()
plotly::save_image(sankeyp, "results/figures/cluster04_DEGs_sankey.png",
  width = 480)
