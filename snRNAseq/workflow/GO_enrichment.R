#################### README
# Aim: run enrichment analyses on veh ccpe DEGs
# Tutorial: https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
# msigdb tutorial: https://yulab-smu.top/biomedical-knowledge-mining-book/027-universal-enrichment.html

#################### Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(gridExtra)
library(org.Hs.eg.db)
library(viridis)
library(ggtree)
library(DOSE)
library(enrichplot)

set.seed(100)

#################### Define variables
{
DEGlist <- "results/DGE/ALNAIK/ALNAIK_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#53d0f7" # CLdn15-high, originally 2, 01 in manuscript
testcol <- "#8361d3" #Stra6-high originally 3, 03 in manuscript

top_up <- c()

}

{
DEGlist <- "results/DGE/ALNASL/ALNASL_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#8361d3" #Stra6-high originally 3, 03 in manuscript
testcol <- "#00ad96" #Stra6-low originally 4, 04 in manuscript

top_up <- c()

}

{
DEGlist <- "results/DGE/SYRMA/SYRMA_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#1fb2df"
testcol <- "#3f6ff5"

top_up <- c("ephrin receptor signaling pathway",
            "cardiac septum morphogenesis",
            "heart valve development",
            "blood vessel development",
            "anatomical structure formation involved in morphogenesis",
            "cell fate commitment",
            "mesenchyme development",
            "regulation of epithelial to mesenchymal transition",
            "collagen-containing extracellular matrix",
            "negative regulation of cell motility",
            "establishment of cell polarity")

}

{
DEGlist <- "results/DGE/SULAFAT/SULAFAT_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#1fb2df"
testcol <- "#8854d0"

# For SULAFAT, determined after running this whole script
top_up <- c("translation",
            "cell activation",
            "cytokine production",
            "cellular response to cytokine stimulus",
            "regulation of tissue remodeling",
            "collagen-containing extracellular matrix",
            "wound healing",
            "actin cytoskeleton organization",
            "regulation of proteolysis",
            "regulation of cell migration")
}

{ # SUBRA: # tst: "17" ref: "05"
DEGlist <- "results/DGE/SUBRA/SUBRA_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#1fb2df"
testcol <- "#255921"
}

{ # SUALOCIN: #tst: "17" ref: "09"
DEGlist <- "results/DGE/SUALOCIN/SUALOCIN_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#3f6ff5"
testcol <- "#255921"
}

{ # STEROPE: #tst: "17", ref: "11"
DEGlist <- "results/DGE/STEROPE/STEROPE_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#8854d0"
testcol <- "#255921"
}

{ # SPICA CLuster 04 (EPDC cluster from main UMAP): Dev vs regen 10dpci
DEGlist <- "results/DGE/SPICA/SPICA04_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#C0EBA3"
testcol <- "#F395B0"

# injury-stuff
top_up <- c("complement activation",
            "MAP kinase phosphatase activity",
            "fatty acid metabolic process",
            "cytokine production",
            "lysosomal membrane"
            )

# devo stuff
top_down <- c("regulation of stem cell proliferation",
              "cytoplasmic translation",
              "regulation of extent of cell growth",
              "transmembrane receptor protein kinase activity",
              "ephrin receptor activity",
              "embryonic pattern specification",
              "cellular response to retinoic acid",
              "cell fate commitment")

}

# 
{ # SITULA CLuster 04 (EPDC cluster from main UMAP): Dev vs regen 14dpci
DEGlist <- "results/DGE/SITULA/SITULA04_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#C0EBA3"
testcol <- "#FFCC66"

# injury stuff
top_up <- c("reactive oxygen species biosynthetic process",
            "MAP kinase phosphatase activity",
            "cytokine production",
            "ligand-gated monoatomic cation channel activity",
            "regulation of synaptic transmission, glutamatergic") #PTGS2/NRXN1/SLC38A2/SERPINE2
# devo stuff
top_down <- c("cytoplasmic translation",
              "cytokine activity",
              "embryonic pattern specification",
              "ephrin receptor activity",
              "mRNA 5'-UTR binding", # LARP1/SYNCRIP/IGF2BP3/RPL5/CCT5/RSL1D1/IGF2BP2/RPS13/GNL3/RPS3A/RPS14/NCL/IGF2BP1
              "mesenchymal cell proliferation",
              "regulation of extent of cell growth",
              "lateral sprouting from an epithelium") #BMP7/FGF10/WNT5A/FGFR2
              # "secondary heart field specification") #WNT5A/WNT11,
              # too few genes maybe? it didn't appear on the ridge plot
}

{ # TEJAT: 03 Dev vs regen 10dpci
DEGlist <- "results/DGE/TEJAT/TEJAT03_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#C0EBA3"
testcol <- "#F395B0"

# injury-stuff
top_up <- c("interleukin-6 production", "peptide hormone binding",
            "G protein-coupled peptide receptor activity",
            "lipid biosynthetic process",
            "inflammatory response")

# devo stuff
top_down <- c("cytoplasmic translation", "chromosome organization",
              "cell cycle DNA replication", "positive regulation of RNA splicing",
              "mRNA 5'-UTR binding", "telomere maintenance",
              "mRNA splicing, via spliceosome")

}

{ # THABIT: 03 Dev vs regen 14dpci
DEGlist <- "results/DGE/THABIT/THABIT03_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#C0EBA3"
testcol <- "#FFCC66"

# injury-stuff
top_up <- c("regulation of fibroblast growth factor receptor signaling pathway",
            "extracellular matrix organization",
            "G protein-coupled peptide receptor activity",
            "interleukin-6 production",
            "response to oxygen levels")

# devo stuff
top_down <- c("cytoplasmic translation",
              "chromosome organization",
              "cell cycle DNA replication",
              "DNA conformation change",
              "mRNA 5'-UTR binding",
              "telomere maintenance")

}


{ # TEJAT: 11 Dev vs regen 10dpci
DEGlist <- "results/DGE/TEJAT/TEJAT11_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#C0EBA3"
testcol <- "#F395B0"

# injury-stuff
top_up <- c("carboxylic acid transport", "amino acid catabolic process",
            "activation of innate immune response", "vacuolar membrane")

# devo stuff
top_down <- c("cytoplasmic translation", "complex of collagen trimers",
              "ubiquitin-protein transferase inhibitor activity",
              "positive regulation of signal transduction by p53 class mediator",
              )

}


{ # EPDC_A Ref: sham, Tst: vehicle.
DEGlist <- "results/DGE/ALGORAB/ALGORABEPDC_A_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#ffacf8" # Sham
testcol <- "#fffbc0" # Vehicle

top_up <- c()

}

{ # EPDC_A Ref: Naive vs. Sham
DEGlist <- "results/DGE/ALKALUROPS/ALKALUROPSEPDC_A_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#b7fff9" # Naive
testcol <- "#ffacf8" # Sham

top_up <- c()

}


{ # EPDC_A Ref: veh, Tst: cCpE
DEGlist <- "results/DGE/ALGOL/ALGOLEPDC_A_DGE_results.xlsx"
log2FC_threshold <- 0
refcol <- "#b2ffab" # Vehicle
testcol <- "#fffbc0" # cCpE


top_up <- c()

}


######################## Specify filenames for results
set.seed(100)

my_palette <- colorRampPalette(plasma(20))(100) # custom viridis colors
options(enrichplot.colours = my_palette) # set viridis as colorscale

path2save <- dirname(DEGlist)
filename <- gsub(pattern = ".xlsx", replacement = "_enrichment.pdf",  x = basename(DEGlist))

filename2 <- gsub(pattern = ".xlsx",
                  replacement = "_enrichment.rds",
                  x = basename(DEGlist))

filename3 <- gsub(pattern = ".xlsx",
                  replacement = "_enrichment.xlsx",
                  x = basename(DEGlist))

filename4 <- gsub(pattern = ".xlsx",
                  replacement = "_GSEA_custom.png",
                  x = basename(DEGlist))

path2save

######################### Reload saved data

# prerun_res <- readRDS(file.path(path2save, filename2))
# GSEA_res <- prerun_res[[5]] # skip to the GSEA section

########################### Load gene symbol mapping file
# Load Zeyu's file
identifiers <- read.csv("../../identifiers/NCBI_aPleWal1.hap1.20221129_genomic_GeneIDtoFeature_20250601_HsandMm_shared.tsv", sep = "\t")
colnames(identifiers)

identifiers <- identifiers %>%
  dplyr::select(c("genename_Hs_PwNameFormat", "genename_Hs"))

############################ Load list of DEGs

# Load DEG list
deglist <- read_xlsx(DEGlist, sheet = "sig_degs")

# Start PDF file to save results
pdf(file.path(path2save, filename))

# Map human symbols to gene identifiers
deglist1 <- left_join(deglist, identifiers,
                      join_by(gene == genename_Hs_PwNameFormat))

########################### Visualize Log2FC value distribution
hist(deglist1$avg_log2FC, 
     main = "log2FC values for significant DEGs",
     col = "darkgreen")

########################### Log2 FC cutoff()
nrow(deglist1) # no log2FC cutoff

deglist2 <- deglist1 %>%
  dplyr::filter(abs(avg_log2FC) > log2FC_threshold)

nrow(deglist2)

########################### Check & deal with duplicate names
# Average the numerical values for Pw genes with same human gene name
# and concatenate the non-numerical columns that have been collapsed
# (like the Pw gene names)

deglist3 <- deglist2 %>%
  group_by(genename_Hs) %>%
  summarise(across(where(is.numeric),
                   mean,
                   na.rm = TRUE),
            across(where(~ !is.numeric(.)),
                   ~ paste(unique(na.omit(.)), collapse = "; ")),
            .groups = "drop")

deglist3 %>% View()

########################## Unranked list of sig & filtered genes
# Extract genenames & log2FC values
filtlist <- deglist3$avg_log2FC
names(filtlist) <- deglist3$genename_Hs

length(filtlist)

########################### Ranked list of all genes tested (including non-sig)
deglist_full <- read_xlsx(DEGlist, sheet = "full_degs")
nrow(deglist_full)

# Map human symbols to gene identifiers
deglist_full1 <- left_join(deglist_full, identifiers,
                      join_by(gene == genename_Hs_PwNameFormat))

# Visualize log2FC distribution
hist(deglist_full1$avg_log2FC,
     main = "log2FC values for all genes tested",
     col = "pink")

# Average the log2FC values
deglist_full2 <- deglist_full1 %>%
  group_by(genename_Hs) %>%
  summarise(across(where(is.numeric),
                   mean,
                   na.rm = TRUE),
            across(where(~ !is.numeric(.)),
                   ~ paste(unique(na.omit(.)), collapse = "; ")),
            .groups = "drop")

deglist_full2 %>% View()

rankedlist <- deglist_full2$avg_log2FC
names(rankedlist) <- deglist_full2$genename_Hs

# sort in decreasing order for ranked analysis
rankedlist <- sort(rankedlist, decreasing = TRUE)

length(rankedlist)
########################### Run overrepresentation analyses
filtlist_up <- filtlist[filtlist > 0]
filtlist_down <- filtlist[filtlist < 0]

length(filtlist_up)
length(filtlist_down)

# unranked all
all_GO <- enrichGO(gene = names(filtlist),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "ALL",
                   minGSSize = 5,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
as.data.frame(all_GO) %>% View()
dotplot(all_GO) + ggtitle("All DEGs GO ORA")

# unranked upregulated
upreg_GO <- enrichGO(gene = names(filtlist_up),
                     OrgDb = org.Hs.eg.db,
                     keyType = "SYMBOL",
                     minGSSize = 5,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05)
as.data.frame(upreg_GO) %>% View()
dotplot(upreg_GO) + ggtitle("Upreg DEGs GO ORA")

# unranked downregulated
downreg_GO <- enrichGO(gene = names(filtlist_down),
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       minGSSize = 5,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05)
as.data.frame(downreg_GO) %>% View()
dotplot(downreg_GO) + ggtitle("Downreg DEGs GO ORA")

########################## Run enrichment analysis
# ranked
GSEA_res <- gseGO(geneList = rankedlist,
                  ont = "ALL",
                  keyType = "SYMBOL",
                  minGSSize = 5,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  OrgDb = org.Hs.eg.db,
                  verbose = FALSE)

as.data.frame(GSEA_res) %>% View()
dotplot(GSEA_res, x = "NES") + ggtitle("All DEGs GSEA")
dotplot(GSEA_res, x = "NES", split = ".sign", showCategory = 15) +
  facet_grid(~.sign)

############################## Visualize GSEA results

# Define terms to visualize: top terms
top_up <- GSEA_res@result %>%
  as.data.frame() %>%
  filter(NES > 0) %>%
  # arrange(NES) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  dplyr::pull(Description)

# bottom terms
top_down <- GSEA_res@result %>%
  as.data.frame() %>%
  filter(NES < 0) %>%
  # arrange(NES) %>%
  arrange(p.adjust) %>%
  slice_tail(n = 10) %>%
  dplyr::pull(Description)

# Color by p-value
options(enrichplot.colours = my_palette) # set viridis as colorscale
# ridgeplot(GSEA_res, showCategory = 175)
ridgeplot(GSEA_res, showCategory = c(top_up, top_down))

# Color by NES value (pos = high in test, neg = high in ref)
options(enrichplot.colours = c(refcol, testcol)) # set viridis as colorscale
# ridgeplot(GSEA_res, showCategory = , fill = "NES")
ridgeplot(GSEA_res, showCategory = c(top_up, top_down), fill = "NES") 
# +
  # coord_cartesian(ylim = c(1.5, 14.6))

ggsave(file.path(path2save, filename4))

############################# Save images
dev.off()

############################# Save RDS

saveRDS(list(deglist3, all_GO, upreg_GO, downreg_GO, GSEA_res),
        file = file.path(path2save, filename2))

############################# Save summary & excel

No_of_genes <- c(nrow(deglist), # number of sig degs
  nrow(deglist1), # after merging with human identifiers
  nrow(deglist2), # after log2FC cutoff, if any
  nrow(deglist3), # after merging genes with same human gene name
  length(filtlist_up), # no. of degs up (i.e. high in test)
  length(filtlist_down), # no. of degs down (i.e. high in ref)
  nrow(deglist_full), # all genes tested
  nrow(deglist_full1), # after merging with human identifiers
  nrow(deglist_full2), # after merging genes with same human gene name
  sum(rankedlist > 0),
  sum(rankedlist <= 0)
)

Comment <- c("A1. Number of sig degs",
  "A2. After merging with human identifiers",
  "A3. After log2FC cutoff, if any",
  "A4. After merging genes with same human gene name",
  "A5. No. of sig degs up (i.e. high in test)",
  "A6. No. of sig degs down (i.e. high in ref)",
  "B1. All genes tested in DGE analyis",
  "B2. After merging with human identifiers",
  "B3. After merging genes with same human gene name",
  "B4. No. of degs up (i.e. high in test)",
  "B5. No. of degs down (i.e. high in ref)"
)

summary_nums <- data.frame(Comment, No_of_genes)

write_xlsx(list(Summary = summary_nums,
                Unranked_sig_genelist = deglist3,
                Ranked_full_genelist = deglist_full2,
                GO_ORA = as.data.frame(all_GO),
                GO_ORA_up = as.data.frame(upreg_GO),
                GO_ORA_down = as.data.frame(downreg_GO),
                GSEA = as.data.frame(GSEA_res)),
           file.path(path2save, filename3))

########################## Get the genes associated with selected GO terms
terms <- c(top_up, top_down)
terms

x <- GSEA_res@result %>%
  dplyr::filter(Description %in% terms) %>%
  dplyr::pull(core_enrichment)

genes <- strsplit(x, "/") %>% unlist()

genes_Pw <- identifiers %>%
  dplyr::filter(genename_Hs %in% genes) %>%
  dplyr::pull(genename_Hs_PwNameFormat)

genes_Pw
