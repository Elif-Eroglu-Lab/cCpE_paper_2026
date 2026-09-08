###################### README
# Author: Inika Prasad, inika.prasad@gmail.com
# Cluster A vs Cluster B, for each condition
# Each column contains log2FC values from a file

##################### Process
# Read in data for each set of DGE results
# file names that start with the unique code and end with results.
# List Excel files
# Make heatmaps for DEGs in a cluster
# veh 
##################### Load libraries
library(readxl)
library(dplyr)

################### OPTION 1: Define the folder and files ######################
{
unique_code <- "ALHENA"
unique_code <- "ALIOTH"
}


{ # Cluster 03 vs Cluster 07, in all conditions (I think)
unique_code <- "ALIFA"
order_conditions <- c("naive", "14dpsham", "10dpci_veh", "10dpci_cCpE")
}

{ # Sham vs. Cryoinjury (veh) for mesothelial epi
unique_code <- "ALGORAB"
order_conditions <- c("03", "07", "04", "13", "14")
}

{ # Sham vs. Cryoinjury, Sham vs cCpE, and veh vs ccpe for meso epi
unique_code <- c("ALGORAB", "ALMAAK", "ALGOL")
order_conditions <- c("ALGORAB03", "ALMAAK03", "ALGOL03",
                      "ALGORAB07", "ALMAAK07", "ALGOL07",
                      "ALGORAB04", "ALMAAK04", "ALGOL04")
# Define annotations
Comparison <- c("Sham vs. Vehicle", "Sham vs. c-CpE", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Sham vs. c-CpE", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Sham vs. c-CpE", "Vehicle vs. c-CpE")
Cluster <- c("Cluster 03", "Cluster 03", "Cluster 03",
             "Cluster 07", "Cluster 07", "Cluster 07",
             "Cluster 04", "Cluster 04", "Cluster 04")

# Define colors for annotations
ann_colors <- list(Comparison = c("Sham vs. Vehicle" = "#646464",
                                  "Sham vs. c-CpE" = "#d6d6d6",
                                  "Vehicle vs. c-CpE" = "#ffb0de"),
                   Cluster = c("Cluster 03" = "#53d0f7",
                               "Cluster 04" = "#8361d3",
                               "Cluster 07" = "#3f6ff5"))
mygaps <- c(3, 6)
}

{ # Sham vs. Cryoinjury, and veh vs ccpe for meso epi
unique_code <- c("ALGORAB", "ALGOL")
order_conditions <- c("ALGORAB03", "ALGOL03",
                      "ALGORAB07", "ALGOL07",
                      "ALGORAB04", "ALGOL04")

# Define annotations
Comparison <- c("Sham vs. Vehicle", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Vehicle vs. c-CpE")
Cluster <- c("Cluster 03", "Cluster 03",
             "Cluster 07", "Cluster 07",
             "Cluster 04", "Cluster 04")

# Define colors for annotations
ann_colors <- list(Comparison = c("Sham vs. Vehicle" = "#AA7974",
                                  "Vehicle vs. c-CpE" = "#AF1D77"),
                   Cluster = c("Cluster 03" = "#6BAED6",
                               "Cluster 04" = "#FC9272",
                               "Cluster 07" = "#6A51A3"))
mygaps <- c(2, 4)
}


{ # Sham vs. Cryoinjury, and veh vs ccpe for epi peelout clusters with DEGs
unique_code <- c("ALGORAB", "ALGOL")
order_conditions <- c("ALGORAB04", "ALGOL04",
                      "ALGORABEPDC_A", "ALGOLEPDC_A",
                      "ALGORAB02", "ALGOL02")

# Define annotations
Comparison <- c("Sham vs. Vehicle", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Vehicle vs. c-CpE",
                "Sham vs. Vehicle", "Vehicle vs. c-CpE")
Cluster <- c("Cluster 04", "Cluster 04",
             "Cluster EPDC_A", "Cluster EPDC_A",
             "Cluster 02", "Cluster 02")

# Define colors for annotations
ann_colors <- list(Comparison = c("Sham vs. Vehicle" = "#dfffa2",
                                  "Vehicle vs. c-CpE" = "#f5c4ff"),
                   Cluster = c("Cluster EPDC_A" = "#1c368d",
                               "Cluster 04" = "#8361d3",
                               "Cluster 02" = "#e75fb3"))
mygaps <- c(2, 4)
}

# USED IN MAIN FIGURE
{ # Sham vs. Cryoinjury, and veh vs ccpe for cluster 04
unique_code <- c("ALGORAB", "ALGOL")
order_conditions <- c("ALGORAB04", "ALGOL04")

# Define annotations
Comparison <- c("Sham vs. Vehicle", "Vehicle vs. c-CpE")
Cluster <- c("Cluster 04", "Cluster 04")

# Define colors for annotations
ann_colors <- list(Comparison = c("Sham vs. Vehicle" = "#dfffa2",
                                  "Vehicle vs. c-CpE" = "#f5c4ff"),
                   Cluster = c("Cluster 04" = "#8361d3"))
mygaps <- c()
}

###################
{ # Megaminus annotations
# Define colors for annotations
ann_colors <- list(Comparison = c("Naive vs. Early Active" = "#FFCC66", # TARAZED
                                  "Sham vs. Early Active" = "#C0EBA3", # TEGMINE
                                  "Early Active vs. 10dpci" = "#94D3D7", # TEJAT
                                  "Early Active vs. 14dpci" = "#F395B0"), # THABIT
                   Cluster = c("Cluster 03" = "#e75fb3",
                               "Cluster 05" = "#1fb2df",
                               "Cluster 09" = "#3f6ff5",
                               "Cluster 11" = "#8854d0",
                               "Cluster 15" = "#5299B8",
                               "Cluster 16" = "#f8d442"))
mygaps <- c()
}

{ # Megaminus: characterizing the Early active (in comp to naive & 14dpsham)
  unique_code <- c("TEGMINE", "TARAZED")
  order_conditions <- c("TARAZED03", "TARAZED05", "TARAZED09", "TARAZED11", "TARAZED15", "TARAZED16",
                        "TEGMINE03", "TEGMINE05", "TEGMINE09", "TEGMINE11", "TEGMINE15", "TEGMINE16")

  # Define annotations
  Comparison <- c(rep("Naive vs. Early Active", 6),
                  rep("Sham vs. Early Active", 6))
  Cluster <- c("Cluster 03", "Cluster 05", "Cluster 09", "Cluster 11", "Cluster 15", "Cluster 16",
              "Cluster 03", "Cluster 05", "Cluster 09", "Cluster 11", "Cluster 15", "Cluster 16")

  log2FC_cutoff <- 2
  filename <- "results/figures_megaminus/TARAZED_TEGMINE.pdf"
  log2FC_dist_filename <- "results/figures_megaminus/TARAZED_TEGMINE_log2FCdistributions.png"
}

#################
{ # Megaminus: characterizing the Injury conditions (in comp to Early active)
  unique_code <- c("TEJAT", "THABIT")
  order_conditions <- c("TEJAT03", "TEJAT05", "TEJAT09", "TEJAT11", "TEJAT15", "TEJAT16",
                        "THABIT03", "THABIT05", "THABIT11", "THABIT15", "THABIT16")

  # Define annotations
  Comparison <- c(rep("Early Active vs. 10dpci", 6),
                  rep("Early Active vs. 14dpci", 6))
  Cluster <- c("Cluster 03", "Cluster 05", "Cluster 09", "Cluster 11", "Cluster 15", "Cluster 16",
              "Cluster 03", "Cluster 05", "Cluster 11", "Cluster 15", "Cluster 16")

  log2FC_cutoff <- 2
  filename <- "results/figures_megaminus/TEJAT_THABIT.pdf"
  log2FC_dist_filename <- "results/figures_megaminus/TEJAT_THABIT_log2FCdistributions.png"
}

#############
{ # Megaminus: characterizing Cluster 03: Injury conditions (in comp to Early active)
  unique_code <- c("TEJAT", "THABIT")
  order_conditions <- c("TEJAT03", "THABIT03")

  # Define annotations
  Comparison <- c(rep("Early Active vs. 10dpci", 1),
                  rep("Early Active vs. 14dpci", 1))
  Cluster <- c("Cluster 03", "Cluster 03")

  log2FC_cutoff <- 2
  filename <- "results/figures_megaminus/Cluster03_TEJAT_THABIT.pdf"
  log2FC_dist_filename <- "results/figures_megaminus/Cluster03_TEJAT_THABIT_log2FCdistributions.png"
}

{ # Megaminus: characterizing Cluster 11: Injury conditions (in comp to Early active)
unique_code <- c("TEJAT", "THABIT")
order_conditions <- c("TEJAT11", "THABIT11")

# Define annotations
Comparison <- c("Early Active vs. 10dpci", "Early Active vs. 14dpci")
Cluster <- c("Cluster 11", "Cluster 11")

log2FC_cutoff <- 2
filename <- "results/figures_megaminus/Cluster11_TEJAT_THABIT.pdf"
log2FC_dist_filename <- "results/figures_megaminus/Cluster11_TEJAT_THABIT_log2FCdistributions.png"
}


{ # Megaminus: characterizing Cluster 11: in Early active (compared to naive & sham)
unique_code <- c("TEGMINE", "TARAZED")
order_conditions <- c("TEGMINE11", "TARAZED11")

# Define annotations
Comparison <- c(rep("Naive vs. Early Active", 1),
                rep("Sham vs. Early Active", 1))
Cluster <- c("Cluster 11", "Cluster 11")

log2FC_cutoff <- 2
filename <- "results/figures_megaminus/Cluster11_TEGMINE_TARAZED.pdf"
log2FC_dist_filename <- "results/figures_megaminus/Cluster11_TEGMINE_TARAZED_log2FCdistributions.png"
}


####################### Get the relevant files
# Make the folder names
folder1 <- paste0("results/DGE/", unique_code)
folder1

# List the files in the folder
files1 <- list.files(folder1, pattern = "\\.xlsx$", full.names = TRUE)
files1

# Remove the summary file
files1 <- files1[!files1 %in% paste0("results/DGE/", unique_code, "/summary_DGE_results.xlsx")]
files1

################## Collect log2FC information in the correct order
# Empty list to collect information about DEGs
all_sig_degs <- list()

# For every condition, collect the significantly changing genes
for (i in seq_along(order_conditions)) {

    # Determine which file has info about that condition
    degname <- order_conditions[i]
    degfile <- files1[grep(degname, files1)]

    # if there actually is a file
    if (length(degfile) != 0) {
        # read in the significant degs for that file/condition
        sig_degs <- read_xlsx(path = degfile,
                              sheet = "sig_degs") %>%
        as.data.frame() %>%
        dplyr::filter(abs(avg_log2FC) > log2FC_cutoff) %>%
        dplyr::select(c(gene, avg_log2FC)) %>%
        # rename every column except the gene column
        rename_with(~ paste0(degname), .cols = -gene)

    all_sig_degs[[degname]] <- sig_degs
  }


}

################## Collect log2FC & adjusted p-val information in the correct order
# Get list of genes that will be used
req_genes <- purrr::reduce(all_sig_degs, full_join, by = "gene") %>%
  dplyr::pull(gene) %>%
  unique()

head(req_genes)
length(req_genes)

# Empty list to collect information about DEGs
all_sig_degs_adjpvals <- list() # to collect adjusted p value matrix
all_sig_degs_fcvals <- list() # to collect DEGs

# For every condition
for (i in seq_along(order_conditions)) {

    # Determine which file has info about that condition
    degname <- order_conditions[i]
    degfile <- files1[grep(degname, files1)]

    # if there actually is a file
    if (length(degfile) != 0) {
        # read in the required degs for that file/condition
        reqdegs_df <- read_xlsx(path = degfile,
                                sheet = "full_degs") %>%
        as.data.frame() %>%
        filter(gene %in% req_genes)
        
        degs_fcvals <- reqdegs_df %>%
        dplyr::select(c(gene, avg_log2FC)) %>%
        # rename every column except the gene column
        rename_with(~ paste0(degname), .cols = -gene)

        degs_adjpvals <- reqdegs_df %>%
        dplyr::select(c(gene, p_val_adj)) %>%
        mutate(adjpvalstars = case_when(p_val_adj < 0.001 ~ "***",
                                p_val_adj < 0.01  ~ "**",
                                p_val_adj < 0.05  ~ "*",
                                TRUE ~ "")) %>%
        dplyr::select(c(gene, adjpvalstars)) %>%
        # rename every column except the gene column
        rename_with(~ paste0(degname), .cols = -gene)

    all_sig_degs_adjpvals[[degname]] <- degs_adjpvals
    all_sig_degs_fcvals[[degname]] <- degs_fcvals
    }


}


all_sig_degs_fcvals <- purrr::reduce(all_sig_degs_fcvals, full_join, by = "gene")
all_sig_degs_adjpvals <- purrr::reduce(all_sig_degs_adjpvals, full_join, by = "gene")
rownames(all_sig_degs_adjpvals) <- all_sig_degs_adjpvals$gene

# QC that we have the correct number of DEGs
dim(all_sig_degs_fcvals)
# dim(all_degs)
dim(all_sig_degs_adjpvals)

######
# Optional filtering: only include genes which changed significantly in veh-ccpe
{ # NEED TO UPDATE
all_sig_degs_adjpvals <- all_sig_degs_adjpvals %>%
  filter(if_any(contains("ALGOL"),
    ~ .x %in% c("*", "**", "***")
  ))

all_sig_degs_fcvals <- all_sig_degs_fcvals[all_sig_degs_fcvals$gene %in% all_sig_degs_fcvals_pvals$gene, ]
}


###################### OPTION 2: Different way to get DEGs #########################

# ALBALI: veh_01 vs. ccpe_01
# ALDEBARAN: veh_06 vs. ccpe_06

# ALBIREO: veh 02 vs ccpe 02
# ALCHIBA: veh_03 vs ccpe_03
# ALCOR: veh_04 vs. ccpe_04
# ALCYONE: veh_05 vs. ccpe_05

################## Specify the DEGs you're interested in
degkeys <- list(Epi1 = "ALBALI",
                Epi2 = "ALDEBARAN",
                EPDC2 = "ALBIREO",
                EPDC3 = "ALCHIBA",
                EPDC4 = "ALCOR",
                EPDC5 = "ALCYONE")

##############
degkeys <- list(ShamActivatedEpi_vs_NaiveEpi1 = "ALKIDR", # tst: "14dpsham_04", ref: "PM_naive_03"
                VehActivatedEpi_vs_NaiveEpi1 = "ALKES") #tst: "10dpci_veh_04", ref: "PM_naive_03"



#################### Read in log2FC and p-values for every DEG
degkeys

all_degs <- list()

for (i in seq_along(degkeys)) {

    degcode <- degkeys[i]
    degname <- names(degkeys[i])

    sig_degs <- read_xlsx(path = paste0("results/DGE/",
                                        degcode, "/",
                                        degcode,
                                        "_DGE_results.xlsx"),
                        sheet = "sig_degs") %>%
    as.data.frame() %>%
    select(c(gene, avg_log2FC)) %>%
    rename_with(~ paste0(degname), .cols = -gene)

    all_degs[[degname]] <- sig_degs
}

all_degs

##################################################################

# Enforce a log2FC cutoff?
all_degs <- all_degs %>%
  mutate(across(where(is.numeric), ~ ifelse(. > -1 & . < 1, 0, .)))
###################### Use all_degs to make a heatmap
# go to a different script for that