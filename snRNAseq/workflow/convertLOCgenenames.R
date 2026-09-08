## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(readxl)
library(ggplot2)

## Read in identifier dataset
identifiers <- read_excel("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.xlsx")

## Read in the Seurat object
sfi <- readRDS(snakemake@input[["rds"]])

# Current genes in Seurat object
current_genes <- rownames(sfi)

# Mapping between identifiers with LOC genes (og_names) and human-readable genes (new_names)
og_identifiers <- identifiers$gene_id_dash
new_identifiers <- identifiers$genename_Hs_PwNameFormat
name_map <- setNames(new_identifiers, og_identifiers)

# Check whether our current genes match the original genes
matched_genes <- intersect(current_genes, og_identifiers)
matched_genes %>% unique() %>% length()
current_genes %>% unique() %>% length()

# Map to the new gene names
new_genenames <- current_genes
new_genenames[match(matched_genes, current_genes)] <- name_map[matched_genes]

# Set the new gene names
rownames(sfi) <- new_genenames

# Save the Seurat object
saveRDS(sfi, snakemake@output[["rds"]])
