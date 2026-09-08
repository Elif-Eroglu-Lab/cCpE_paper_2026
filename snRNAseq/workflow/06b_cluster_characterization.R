# Setup log file
logfile <- file(snakemake@log[[1]], open = "wt")
sink(logfile)                          # Redirect stdout
sink(logfile, type = "message")       # Redirect stderr

# Aim: run enrichment analyses on cluster markers
# Suitable for V2 genome aligned samples as it uses the Trinotate custom GO database

## Load libraries
library(Seurat)
library(dplyr)
library(writexl)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(gridExtra)

set.seed(100)

## Define snakemake variables
markerlist <- snakemake@input[["markerlist"]]
dir4results <- snakemake@output[["dotplot"]] %>% dirname()
filename <- snakemake@output[["dotplot"]] %>% basename()
whichont <- "all"
ontologyfile <- "miscellaneous/GO_all.csv"
dotplot_combo <- snakemake@output[["dotplot_combo"]]
tester <- snakemake@output[["tester"]]

ncol <- 2
nrow <- 3

print(markerlist)
print(dir4results)
print(filename)
print(ontologyfile)
print(dotplot_combo)
ont <- read.csv(ontologyfile)

## Function to run enrichment analyses on a list of gene lists
#' @title newtpathwayORA_lists
#' @description do pathway overrepresentation analysis using enrichGO, save results conveniently.
#' @param listof_featurelists A list containing lists of genes. Each list will go through gene enrichment. Each list must have a name, because it will be used as a label
#' @param ontology select ontology for the GO enrichment. options are "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular compartment, and "all"
#' @param savedotplot boolean for whether to save the resulting dotplots
#' @param saveexcel boolean for whether to save enrichment results in an excel file  
#' @param dir2save directory in which results should be saved
#' @param list_titles a list of what each genelist is. The dotplot
newtpathwayORA_lists <- function(listof_featurelists,
                                 ontology = "all",
                                 savedotplot = TRUE,
                                 dotplot_combo = "combodotplot.pdf",
                                 saveexcel = TRUE,
                                 dir2save = NULL,
                                 filename = "newtGOdotplot.pdf",
                                 nrow = 3,
                                 ncol = 3,
                                 ontology_file = "miscellaneous/Pw_genome2/GO_all.csv"){

  if (is.null(listof_featurelists)) (stop("No features provided"))
  if (is.null(names(listof_featurelists))) (stop("Feature lists must be named"))

  #Load data for appropriate ontology
  if (!file.exists(ontology_file)) (stop("Ontology file does not exist in specified location"))
  ont <- read.csv(ontology_file)

  # Save the dotplots if specified
  if (savedotplot == TRUE) {
    dotplot_path <- file.path(dir2save, filename)
    plot_list <- list()
    pdf(dotplot_path)
  }

  # Loop through the provided lists
  for (i in seq_along(listof_featurelists)){

    label <- names(listof_featurelists)[i]
    featurelist <- listof_featurelists[[i]] %>% unlist() %>% as.vector()
    print(label)

    # Run ORA on every set of genes in the lists
    ORA_results <- enricher(gene = featurelist,
                            TERM2GENE = ont,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            universe = NULL)

    # TO DO : skip dotplot if there aren't enough significant terms

    # Make dotplot
    p1 <- dotplot(ORA_results) +
      ggtitle(paste0(label, ": enriched GO ", ontology))
    print(p1)

    p <- p1 +
      theme_minimal(base_size = 8) +
      theme(plot.margin = margin(3, 3, 3, 3))
    plot_list[[i]] <- p


    # Save excel file if specified
    if (saveexcel == TRUE) {
      excel_path <- paste0(dir2save, "/", label, "_newtGO.xlsx")
      write_xlsx(
                 list("ORA_results" = ORA_results@result,
                      "input_genes" = as.data.frame(featurelist)),
                 path = excel_path)
    }
  }
  if (savedotplot == TRUE) {
    dev.off()
    grob_pages <- marrangeGrob(grobs = plot_list, nrow = nrow, ncol = ncol)
    ggsave(dotplot_combo, grob_pages, width = 8.27, height = 11.69, units = "in")
  }
}

# Load marker list
m <- read_xlsx(markerlist)

# Split up the marker list by cluster
cluster_markerlists <- m %>%
  group_by(cluster) %>%
  dplyr::select(gene) %>%
  group_split(.keep = FALSE) %>%
  #name the lists
  setNames(paste0("cluster", m %>%
                    group_by(cluster) %>%
                    group_keys() %>%
                    pull(cluster)))

# Run enrichment on each cluster: individually
newtpathwayORA_lists(listof_featurelists = cluster_markerlists,
                     ontology = whichont,
                     dotplot_combo = dotplot_combo,
                     savedotplot = TRUE,
                     saveexcel = TRUE,
                     dir2save = dir4results,
                     filename = filename,
                     ontology_file = ontologyfile,
                     nrow = nrow,
                     ncol = ncol)

# Close log file
sink(type = "message")
sink()
close(logfile)