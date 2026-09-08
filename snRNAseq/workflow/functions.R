#' @title Plot elbow plot
#' @description Creates elbow plot for a dimensionality reduction
#' @param object The Seurat object for which you want the elbow plot
#' @param dimred The dimensionality reduction. E.g. "pca"
plot_elbow <- function(object, dimred = "pca") {
  pct <- object[[dimred]]@stdev / sum(object[[dimred]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90%
  # and % variation associated with the PC as less than 5
  # to identify when getting more PCs has diminishing returns.
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  # Minimum of the two calculation
  pcs <- min(co1, co2)
  message(paste0("elbow: ", pcs))

  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

  # Elbow plot to visualize
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "grey") +
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_minimal()+
    theme(plot.background = element_rect(fill = "white", color = NA))

  return(p)
}

#' @description Ranks gene from a FindMarkers() output
#' @param m Dataframe with columns for adjusted p value and log2 fold change
#' @param pv Character denoting column name for adjusted p val
#' @param lf Character denoting column name for log 2 fold change
#' @return Returns an ordered dataframe with column score
#'
ranker <- function(m, pv = "p_val_adj", lf = "avg_log2FC") {
    m %>%
        mutate(
          abundance_diff = pct.1 - pct.2, # tst pct - ref pct
          score = -log10(!!sym(pv)) * abs(!!sym(lf)),
          score2 = score * abs(abundance_diff),
          score3 = abs(!!sym(lf)) * abs(abundance_diff),
          score4 =  score * abundance_diff,
          score5 = abs(!!sym(lf)) * abundance_diff) %>%
        arrange(-score) %>%
        return()
}

#' @description Creates a violin plot
#' @param dfr A data.frame with group, gene and value
#' @param ncol Number of columns
#' @param legend_nrow Number of rows of legend
#' @param cols Colours
#' @param base_size Base text size
#'
plot_violin <- function(dfr, ncol = NULL, legend_nrow = NULL, cols = c("#bebada", "#8dd3c7"), base_size = 12) {
  if (is.null(ncol)) ncol <- length(base::unique(dfr$gene))
  ggplot(dfr, aes(x = group, y = value)) +
    geom_violin(aes(colour = group, fill = group, group = group), size = 0.2, scale = "width", width = 1, position = position_dodge(width = 0.8)) +
    facet_wrap(~gene, scales = "free", ncol = ncol, dir = "v") +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    labs(y = "Log normalised expression") +
    guides(fill = guide_legend(nrow = legend_nrow)) +
    coord_flip() +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "top",
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.ticks.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
}


#' @description Creates pval histogram from a FindMarkers() output
#' @param m Dataframe with column p_val
#' @param vp Variable name for p value column
# http://varianceexplained.org/statistics/interpreting-pvalue-histogram/ 
plot_pvaldist <- function(m,vp="p_val") {
  ggplot(m, aes(x = !!sym(vp))) +
    geom_histogram() +
    labs(x = "Raw p value", y = "Counts") +
    theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

#' @description Creates qqplot from a FindMarkers() output
#' @param m Dataframe with columns p_val and p_val_adj
#' @param type Type of p value to plot. "raw", "adjusted" or "both".
#' @param ci Confidence interval
#' @param vp Variable name for p value column
#' @param vpa Variable name for adjusted p value
#' https://github.com/THERMOSTATS/RVA/blob/master/R/plot_qq.R
#' @importFrom ggpointdensity geom_pointdensity
#' @importFrom ggnewscale new_scale_color new_scale_fill
#'
plot_qq <- function(m, type = "both", ci = 0.95, vp = "p_val", vpa = "p_val_adj") {
  
  if (!type %in% c("raw", "adjusted", "both")) stop("Argument 'type' must be one of 'raw', 'adjusted' or 'both'.")
  m <- drop_na(m,!!sym(vp))
  if (nrow(m) == 0) stop("No genes after removing NAs.")
  
  cols_adj <- c("#a3dac9", "#93d4c0", "#84ceb7", "#75c8ae", "#66c2a5", "#5bae94", "#519b84", "#478773", "#3d7463")
  cols_raw <- c("#fcb2aa", "#fca69c", "#fb998e", "#fb8c80", "#fb8072", "#e17366", "#c8665b", "#af594f", "#964c44")
  cols <- c("#66c2a5", "#fb8072")
  names(cols) <- c("adjusted", "raw")
  types <- NULL
  
  p <- ggplot() +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = 3, linewidth = 0.3, color = "grey70")
  
  if (type == "raw" || type == "both") {
    types <- c(types, "raw")
    
    ps <- dplyr::pull(m,!!sym(vp))
    n <- length(ps)
    dfr_raw <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(stats::ppoints(n)),
      clower = -log10(stats::qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper = -log10(stats::qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)),
      type = "raw"
    )
    
    p <- p +
      geom_line(data = dfr_raw, aes(x = expected, y = observed), stat = "smooth", method = "loess", color = "#fb8072", linewidth = 0.4, alpha = 0.2) +
      geom_point(data = dfr_raw, aes(x = expected, y = observed), color = "grey20", fill = "grey20", shape = 20, alpha = 1, size = 1.7, show.legend = FALSE) +
      # geom_point(aes(x=expected,y=observed,color=type,fill=type),shape=20,alpha=1,size=1.3)+
      geom_pointdensity(data = dfr_raw, aes(x = expected, y = observed), shape = 20, alpha = 1, size = 1.5, show.legend = FALSE) +
      # geom_line(data=dfr_raw,aes(expected,cupper),color="#fb8072",size=0.4,linetype=1)+
      geom_ribbon(data = dfr_raw, aes(x = expected, ymin = clower, ymax = cupper, fill = type), alpha = 0.2) +
      scale_color_gradientn(colours = cols_raw)
  }
  
  if (type == "adjusted" || type == "both") {
    types <- c(types, "adjusted")
    if (type == "both") p <- p + new_scale_color()
    
    ps <- dplyr::pull(m,!!sym(vpa))
    n <- length(ps)
    dfr_adj <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(stats::ppoints(n)),
      clower = -log10(stats::qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper = -log10(stats::qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)),
      type = "adjusted"
    )
    
    p <- p +
      geom_line(data = dfr_adj, aes(x = expected, y = observed), stat = "smooth", method = "loess", color = "#66c2a5", size = 0.4, alpha = 0.2) +
      geom_point(data = dfr_adj, aes(x = expected, y = observed), color = "grey20", fill = "grey20", shape = 20, alpha = 1, size = 1.7, show.legend = FALSE) +
      # geom_point(aes(x=expected,y=observed,color=type,fill=type),shape=20,alpha=1,size=1.3)+
      geom_pointdensity(data = dfr_adj, aes(x = expected, y = observed), shape = 20, alpha = 1, size = 1.5, show.legend = FALSE) +
      # geom_line(data=dfr_adj,aes(expected,cupper),color="#66c2a5",size=0.4,linetype=1)+
      geom_ribbon(data = dfr_adj, aes(x = expected, ymin = clower, ymax = cupper, fill = type), alpha = 0.2) +
      scale_color_gradientn(colours = cols_adj)
  }
  
  p <- p +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.4, linetype = 2, color = "black", alpha = 0.5) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.9))) +
    scale_fill_manual(values = cols, drop = TRUE, limits = sort(unique(types))) +
    labs(
      x = expression(paste("Expected -log"[10], " P value")),
      y = expression(paste("Observed -log"[10], " P value"))
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.title = element_blank()) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  return(p)
}

plot_volcano <- function(m,
                         ref_grp,
                         test_grp,
                         vfc = "avg_log2FC",
                         vpa = "p_val_adj",
                         threshold = 0.05,
                         genes2label = NULL,
                         label_size = 3,
                         cols = c("#2166ac", "#b2182b", "#a2a2a2")) {
  
  # Include shortened gene names for ease of reading
  m <- m %>%
    mutate(shortened_gene = sub("-.*", "", gene))

  # Define x axis title
  if (missing(test_grp) | missing(ref_grp)) {
    xlab <- expression(paste0("Mean Log"[2], " fold change"))
  } else {
    xlab <- paste("Mean Log2 fold change\n< Lower  (in ",
                  test_grp, " relative to ", ref_grp, ")  Higher >",
                  sep = "")
  }

  # Add information about direction and significance for genes
  m1 <- m %>%
    mutate(y = -log10(!!sym(vpa)),
           dir = ifelse(!!sym(vfc) > 0 & !!sym(vpa) < threshold, "up",
                        ifelse(!!sym(vfc) < 0 & !!sym(vpa) < threshold, "down",
                               "non-sig")),
           sig = factor(ifelse(!!sym(vpa) < threshold, "sig", "xsig"))) %>%
    mutate(dir = factor(dir, levels = c("up", "down", "non-sig")))

  m2 <- m1

  # Extract genes to label
  if (is.null(genes2label)){
    lg <- pull(slice_head(ranker(m), n = 30), gene)
  } else {
    lg <- genes2label
  }

  m3 <- filter(m2, gene %in% lg)

  # Create volcano plot
  p <- ggplot() +
    geom_point(data = m1,
               aes(x = !!sym(vfc), y = y),
               color = "grey20", fill = "grey20", shape = 21,
               alpha = 1, size = 1.5, show.legend = FALSE) +
    geom_hline(aes(yintercept = -log10(threshold)),
               linetype = 3, linewidth = 0.3, color = "grey70") +
    geom_vline(aes(xintercept = 0),
               linetype = 3, linewidth = 0.3, color = "grey70") +
    geom_point(data = m2,
               aes(x = !!sym(vfc), y = y, fill = dir, color = dir),
               shape = 21, alpha = 1, size = 1.3, show.legend = FALSE) +
    # Add labels
    ggrepel::geom_text_repel(data = m3,
                             aes(x = !!sym(vfc),
                                 y = y,
                                 label = shortened_gene),
                             color = "black",
                             min.segment.length = 0.2,
                             force_pull = 0.8,
                             size = label_size,
                             alpha = 1,
                             segment.size = 0.14,
                             fontface = "bold",
                             show.legend = FALSE) +
    scale_fill_manual(values = setNames(cols,
                                        c("down", "up", "non-sig"))) +
    scale_colour_manual(values = setNames(cols,
                                          c("down", "up", "non-sig"))) +
    labs(x = xlab, y = expression(paste("- Log"[10], " Adjusted p-value"))) +
    theme_classic(base_size = 10) +
    theme(plot.background = element_rect(fill = "white",
                                         color = NA))
  return(p)
}


#' @description Creates violin plot from a seurat object
#' @param obj A seurat object
#' @param vari A variable column name
#' @param cond Levels under the vari variable ordered as ref,test
#' @param g One or more gene names
#' @param assay Assay
#' @param slot Slot
#' @param ncol Number of columns
#' @param legend_nrow Number of rows for legend
#'
plot_violin_seurat <- function(obj,
                               vari, cond = NULL, g, assay = "RNA",
                               slot = "data", ncol = NULL, legend_nrow = NULL,
                               cols = c("#bebada", "#8dd3c7"),
                               ...) {
  if (length(vari) > 1) stop("Argument 'vari' must be of length one.")

  x <- FetchData(obj, c(vari, g), assay = assay, slot = slot)
  if (!is.null(cond)) x <- x %>% dplyr::filter(!!sym(vari) %in% cond)
  x_long <- x %>%
    tidyr::pivot_longer(cols = all_of(g), names_to = "gene", values_to = "value") %>%
    # mutate(value=2^value) %>%
    dplyr::mutate(gene = factor(as.character(gene), levels = g)) %>%
    dplyr::rename(group = !!sym(vari)) %>%
    dplyr::mutate(group = factor(as.character(group), levels = cond))
  plot_violin(
    x_long,
    ncol = ncol,
    legend_nrow = legend_nrow,
    cols = cols,
    ...
  )
}

###############################################################################

#' @title dge
#' @description Runs pairwise DGE
#' @param obj Seurat object
#' @param dfr Comparisons dataframe
#' @param ds Downsample? NULL or a number. Number of cells per identity group.
#' @param export Export?
#' @param path Path for export
#' @param test.use wilcox, MAST etc
#' @param assay Assay
#' @param latent.vars Latent variables
#' @param pvaldist P value distribution histogram
#' @param qqplot P value qqplot.
#' @param volcano Logical exports volcano plot
#' @param violin Logical exports violin plot
#'

dge <- function(obj = NULL, dfr = NULL, ds = 200, export = TRUE, path = NULL, test.use = "MAST", assay = "RNA", latent.vars = NULL, pvaldist = FALSE, qqplot = FALSE, volcano = FALSE, violin = FALSE) {
  if (is.null(obj)) stop("obj missing.")
  if (is.null(dfr)) stop("dfr missing.")
  if (export) {
    if (is.null(path)) stop("path missing.")
  }

  set.seed(100)

  dlist <- vector("list", length = nrow(dfr))
  for (i in 1:nrow(dfr)) {
    cat(paste0(round((i / nrow(dfr)) * 100, 2), "%\n"))

    ref <- dfr$ref[i]
    tst <- dfr$test[i]
    grp <- dfr$grp[i]
    if (!is.na(dfr$covariates[i])){
      latent.vars <- dfr$covariates[i] %>% strsplit(split = ",") %>% unlist() %>%trimws()
    }
    label <- paste0("DGE-", grp, "-", ref, "_vs_", tst, "_", dfr$covariates[i])

    # check if output already exists
    if (!dir.exists(file.path(path, label))) {
      message(paste0("Running grp: ", grp, ", ref: ", ref, ", test: ", tst))
      dir.create(file.path(path, label))
      obj_dge <- obj
      Idents(obj_dge) <- obj_dge[[]][, grp]

      # check if ref and test levels are present in grp
      if (any(obj_dge[[]][, grp] %in% ref) && any(obj_dge[[]][, grp] %in% tst)) {

        # subset
        obj_dge <- subset(obj_dge, subset = !!sym(grp) %in% c(ref, tst))

        # downsampling
        if (!is.null(ds)) {
          if (any(as.data.frame(table(Idents(obj_dge)))$Freq < ds)) {
            message(paste0("Downsampling skipped because one or more groups have cells less than ", ds, " cells."))
            print(setNames(as.data.frame(table(Idents(obj_dge))), c("group", "cells")))
          } else {
            message(paste0("Downsampling data to ", ds, " cells per group."))
            obj_dge <- subset(obj_dge, downsample = ds)
          }
        }

        # get counts for both groups
        nc <- obj_dge[[]] %>%
          dplyr::filter(!!sym(grp) %in% c(ref, tst)) %>%
          dplyr::group_by(!!sym(grp)) %>%
          dplyr::summarise(n = n()) %>%
          dplyr::ungroup() %>%
          tidyr::pivot_wider(names_from = grp, values_from = "n") %>%
          dplyr::select(!!ref, !!tst)

        # if atleast 3 cells
        if (all(unlist(nc) > 2)) {
          m <- FindMarkers(obj_dge,
                           ident.1 = tst, ident.2 = ref, group.by = grp,
                           test.use = test.use, only.pos = FALSE, min.pct = 0.25,
                           logfc.threshold = 0, min.cells.group = 3, random.seed = 100,
                           assay = assay, slot = "data",
                           latent.vars = latent.vars
          ) %>%
            tibble::rownames_to_column("gene") %>%
            dplyr::select(gene, everything()) %>%
            ranker()

          if (export) {
            lb <- file.path(path,label, paste0(label, ".xlsx"))
            writexl::write_xlsx(m, lb, col_names = T)
            message(paste0(lb, " exported."))
          }

          # summary
          dlist[[i]] <- data.frame(
            grp = grp, ref = ref, test = tst,
            cells_ref = unlist(nc)[1], cells_test = unlist(nc)[2],
            degs_up = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC > 0)),
            degs_down = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC < 0)),
            degs_tot = nrow(dplyr::filter(m, p_val_adj < 0.05)),
            skewness = moments::skewness(m$p_val),
            kurtosis = moments::kurtosis(m$p_val),
            label = label,
            statistical_test = test.use, 
            downsampling = ds,
            covariates = paste(latent.vars, collapse = ", "),
            run_on = Sys.time(),
            comment = NA
          )

          if (pvaldist && export) {
            lb <- file.path(path, label, paste0(label, "-pvaldist.png"))
            ggsave(lb, plot_pvaldist(m), height = 4.5, width = 4.5)
            message(paste0(lb," exported."))
          }

          if (qqplot && export) {
            lb <- file.path(path, label, paste0(label, "-qq.png"))
            ggsave(lb, plot_qq(m), height = 4.5, width = 4.5)
            message(paste0(lb," exported."))
          }

          if (volcano && export) {
            lb <- file.path(path, label, paste0(label, "-volcano.png"))
            ggsave(lb, plot_volcano(m, ref_grp = ref, test_grp = tst), height = 4.5, width = 4.5)
            message(paste0(lb, " exported."))
          }

          if (violin && export) {
            m1 <- m %>%
              ranker() %>%
              slice_head(n = 10)
            if (nrow(m1) > 0) {
              p <- plot_violin_seurat(obj_dge, grp, c(ref, tst), g = m1$gene, assay = assay, ncol = 2)
              lb <- file.path(path, label, paste0(label, "-violin.png"))
              ggsave(lb, p, width = 5, height = 7)
              message(paste0(lb, " exported."))
            } else {
              message("Violin plot not exported.")
            }
          }
        } else {
          message("DGE not run. Less than 3 cells.")
        }
      } else {
        message(paste0("One of the levels ", paste(ref, tst, collapse = " or "), " not found in ", grp, "."))
      }
    } else {
      message(paste0("Skipping grp: ", grp, ", ref: ", ref, ", test: ", tst,". File already exists."))
    }
  }
  d <- bind_rows(dlist)
  # rownames(d) <- 1:nrow(d)
  if (export) {
    xlsx::write.xlsx(x = d, file = file.path(path, paste0("dge-summary.xlsx")),  
                     sheetName = Sys.Date, append = TRUE, row.names = FALSE)
    message("dge-summary.xlsx exported.")
  }

  return(d)
}

###############
# Aim: to create a stacked barplot
# group_var = will comprise the colors in each bar. string
# cluster_var = will be on the x axis. string

library(dplyr)
library(ggplot2)
library(rlang)

plot_cluster_proportions <- function(seurat_obj, group_var, cluster_var, sample_name = NULL) {
  # group_var and cluster_var are expected as strings
  df <- seurat_obj@meta.data %>%
    count(!!sym(group_var), Cluster = !!sym(cluster_var)) %>%
    group_by(Cluster) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  
  p <- ggplot(df, aes(x = Cluster, y = pct, fill = !!sym(group_var))) +
    geom_col(color = "black") +
    ylab("Percentage") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    ggtitle(
      label = sample_name %||% "Sample",
      subtitle = paste("Proportions of", group_var, "per cluster")
    )
  
  return(p)
}

####################### Do DGE analysis
dge_function <- function(seurat_object,
                         unique_code,
                         reduction,
                         tst,
                         ref,
                         grp,
                         ds,
                         covariates,
                         test_use,
                         assay,
                         table1,
                         plot1,
                         plot2) {

  # Check if grouping exists in Seurat object
  if (!grp %in% colnames(obj2@meta.data)) {
    stop(paste0(
      grp,
      " is not present in the Seurat object.\n",
      "Available metadata/clusterings are: ",
      paste(colnames(obj2@meta.data), collapse = ", ")
    ))
  }

  available_levels <- obj2@meta.data[[grp]] %>% unique()
  levels2check <- c(ref, tst)

  # Check if ref & tst are present in that grouping
  if (!all(levels2check %in% available_levels)) {
    missing_levels <- levels2check[!levels2check %in% available_levels]
    stop(paste0(
      "The following ref/tst is not present in the metadata column called ",
      grp,
      ": ",
      paste(missing_levels, collapse = ", "), "\n",
      "Available ref/tst are: ",
      paste(available_levels, collapse = ", ")
    ))
  }

  # Parse covariates provided
  if (!is.null(covariates)) {
    latent.vars <- covariates %>% strsplit(split = ",") %>% unlist() %>% trimws()

      # Check if covariates are present in Seurat object
    if (!all(latent.vars %in% colnames(obj2@meta.data))) {
      missing_vars <- latent.vars[!latent.vars %in% colnames(obj2@meta.data)]
      stop(paste0(
        "The following covariates are not present in the metadata column:",
        paste(missing_vars, collapse = ", "), "\n",
        "Available metadata/clusterings are: ",
        paste(colnames(obj2@meta.data), collapse = ", ")
      ))
    }
  } else {
    latent.vars <- NULL
  }


  ################################ DGE from scratch

  description <- paste0("DGE-", grp, "-", ref, "_vs_", tst, "_", covariates)
  print(description)

  # Copy of Seurat object to subset & use for DGE
  Idents(obj2) <- grp
  obj_dge <- obj2

  # Assign identities to Seurat object based on grouping
  Idents(obj_dge) <- grp

  # check if ref and test levels are present in grp
  Idents(obj_dge) %>% unique()

  # Subset object to include test and reference cells
  obj_dge <- subset(obj_dge, subset = !!sym(grp) %in% c(ref, tst))

  # Save this non-downsampled original object for subsequent heatmap generation
  obj_nods <- obj_dge

  # Find number to downsample to, if downsampling 200 is not possible
  ncells_per_grp <- Idents(obj_dge) %>% table %>% as.data.frame()
  colnames(ncells_per_grp) <- c("group", "cells")
  ncells_per_grp


  if (any(ncells_per_grp$cells < ds)) {
    new_ds <- min(ncells_per_grp$cells)
    message(paste0("One or more groups have cells less than ",
                  ds,
                  " cells. Downsampling to the minimum no. of cells instead, 
                  which is ",
                  new_ds))
    print(ncells_per_grp)

    ds <- new_ds

    # Do the downsampling
    message(paste0("Downsampling data to ", ds, " cells per group."))
    obj_dge <- subset(obj_dge, downsample = ds)

  } else if(any(ncells_per_grp$cells >= ds)) {
    # Do the downsampling
    message(paste0("Downsampling data to ", ds, " cells per group."))
    obj_dge <- subset(obj_dge, downsample = ds)

  } else if (is.null(ds)) {
    ds <- "no downsampling"
  }


  # Get count data for both groups (after possible downsampling)
  ncells_used <- Idents(obj_dge) %>% table %>% as.data.frame()
  colnames(ncells_used) <- c("group", "cells")
  ncells_used

  # Must have at least 3 cells per group
  if (sum(ncells_used$cells < 3) != 0) {
    stop("There are less than 3 cells in at least one group")
  }

  # Do the DGE testing
  m <- FindMarkers(obj_dge,
                  ident.1 = tst,
                  ident.2 = ref,
                  group.by = grp,
                  test.use = test_use,
                  only.pos = FALSE,
                  min.pct = 0.25,
                  logfc.threshold = 0,
                  min.cells.group = 3,
                  random.seed = 100,
                  assay = assay,
                  slot = "data",
                  latent.vars = latent.vars) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, everything()) %>%
    ranker()

  run_on <- Sys.time()

  # Subset DEGs to include only significant adj.p.val < 0.05
  sig_degs <- m %>% filter(p_val_adj < 0.05)

  # Construct summary
  dlist <- data.frame(
    unique_code = unique_code,
    grp = grp, ref = ref, test = tst,
    cells_ref = ncells_used[ncells_used$group == ref, ]$cells,
    cells_test = ncells_used[ncells_used$group == tst, ]$cells,
    degs_up = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC > 0)),
    degs_down = nrow(dplyr::filter(m, p_val_adj < 0.05 & avg_log2FC < 0)),
    degs_tot = nrow(dplyr::filter(m, p_val_adj < 0.05)),
    skewness = moments::skewness(m$p_val),
    kurtosis = moments::kurtosis(m$p_val),
    description = description,
    statistical_test = test_use,
    downsampling = ds,
    covariates = paste(latent.vars, collapse = ", "),
    run_on = run_on,
    comment = NA
  )

  # Save info
  writexl::write_xlsx(list(sig_degs = sig_degs,
                           full_degs = m,
                           run_details = dlist),
                      table1)

  # Top 10 genes for violin plots
  m1 <- m %>%
    ranker() %>%
    slice_head(n = 10)

  # Cells used for the DGE in each group
  ref_cells <- colnames(obj_dge)[obj_dge[[grp]] == ref]
  tst_cells <- colnames(obj_dge)[obj_dge[[grp]] == tst]

  # Save result plots in a pdf
  pdf(plot1, height = 7, width = 7)
  p0 <- plot_volcano(m, ref_grp = ref, test_grp = tst)

  p1 <- DimPlot(obj2, cells.highlight = ref_cells, reduction = reduction) +
    ggtitle(ref) +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5))

  p2 <- DimPlot(obj2, cells.highlight = tst_cells, reduction = reduction) +
    ggtitle(tst) +
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5))

  print(p0)
  print(p1)
  print(p2)

  p3 <- plot_pvaldist(m) + ggtitle("Raw p-value distribution")
  p4 <- plot_qq(m) + ggtitle("QQ plot")
  p5 <- plot_violin_seurat(obj_dge, grp,
                           c(ref, tst),
                           g = m1$gene,
                           assay = assay,
                           ncol = 2,
                           base_size = 6)
  print(p3)
  print(p4)
  print(p5)
  dev.off()

  # Make a dotplot showing where the genes are expressed
  if (nrow(sig_degs) > 0) {
    pdf(plot2,
        height = 2.5 + (0.2 * nrow(sig_degs)),
        width = 4 + 0.4 * length(unique(Idents(obj2))))

    p6 <- DotPlot(obj2, features = sig_degs$gene,
                  cols = "RdBu", scale = TRUE) +
      RotatedAxis() +
      coord_flip()

    p7 <- DotPlot(obj2, features = sig_degs$gene,
                  cols = "RdBu", scale = FALSE) +
      RotatedAxis() +
      coord_flip()

    print(p6)
    print(p7)

    dev.off()

  }


}


order_dotplot_genes <- function(dp) {

  # dotplot is a dotplot object from Seurat
  # from which data about percentage expression and average scaled expression
  # can be extracted

  # Percentage expression of genes
  pct_exp_mat <- dp$data %>%
    dplyr::select(features.plot, id, pct.exp) %>%
    pivot_wider(names_from = c(id),
                values_from = c(pct.exp)) %>%
    as.data.frame %>%
    tibble::column_to_rownames("features.plot") %>%
    as.matrix() %>%
    t() %>% scale %>% t()

  # Scaled expression of genes
  scaled_avg_exp_mat <- dp$data %>%
    dplyr::select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(names_from = c(id),
                values_from = c(avg.exp.scaled)) %>%
    as.data.frame %>%
    tibble::column_to_rownames("features.plot") %>%
    as.matrix() %>%
    t() %>% scale %>% t()

  # Combine matrices for pct expression and scaled avg expression
  bigmat <- cbind(pct_exp_mat, scaled_avg_exp_mat)

  # Compute distances between the genes
  gene_dist <- bigmat %>% dist()

  # Cluster genes based on distances
  gene_clust <- hclust(gene_dist, method = "ward.D2")

  # Extract the gene names in order
  ordered_genes <- bigmat[gene_clust$order, ] %>% rownames()

  return(ordered_genes)

}

##########################
barcodeplot_mod <- function (statistics, index = NULL, index2 = NULL, gene.weights = NULL, 
    weights.label = "Weight", labels = c("Down", "Up"), quantiles = c(-1, 
        1) * sqrt(2), col.bars = NULL, alpha = 0.4, worm.max = NULL, worm = TRUE, 
    span.worm = 0.45, xlab = "Statistic", rectangle1col = "lightblue", rectangle2col = "lightpink",...) {
    if (!is.vector(statistics, mode = "numeric")) 
        stop("statistics should be a numeric vector")
    nstat <- length(statistics)
    if (is.null(index)) {
        if (is.null(gene.weights)) {
            stop("Must specify at least one of index or gene.weights")
        }
        else {
            if (length(gene.weights) == nstat) {
                index <- rep_len(TRUE, nstat)
                index2 <- NULL
            }
            else {
                stop("No index and length(gene.weights) doesn't equal length(statistics)")
            }
        }
    }
    else {
        if (anyNA(index)) 
            stop("Need to provide index without NAs")
        if (is.logical(index)) 
            if (length(index) != nstat) 
                stop("Length of index disagrees with statistics")
        if (length(index) > nstat) 
            stop("Length of index disagrees with statistics")
    }
    if (!is.null(index2)) {
        if (!is.null(gene.weights)) 
            warning("gene.weights ignored")
        gene.weights <- statistics
        gene.weights[] <- 0
        gene.weights[index] <- 1
        gene.weights[index2] <- -1
        index <- rep_len(TRUE, nstat)
        index2 <- NULL
    }
    if (!is.null(gene.weights)) {
        if (!is.vector(gene.weights, mode = "numeric")) 
            stop("gene.weights should be a numeric vector")
        if (anyNA(gene.weights)) 
            stop("Need to provide gene.weights without NAs")
        if (all(gene.weights == 0)) 
            stop("gene.weights equal to zero: no selected genes to plot")
        if (length(gene.weights) != length(statistics[index])) 
            stop("Length of gene.weights disagrees with size of set")
        one <- all(gene.weights >= 0) | all(gene.weights <= 0)
        if (one) {
            index2 <- NULL
            gene.weights1 <- rep_len(0, nstat)
            names(gene.weights1) <- names(statistics)
            gene.weights1[index] <- gene.weights
            index <- rep_len(FALSE, nstat)
            names(index) <- names(statistics)
            index[gene.weights1 != 0] <- TRUE
            gene.weights1 <- gene.weights1[gene.weights1 != 0]
            gene.weights <- gene.weights1
        }
        else {
            gene.weights12 <- rep_len(0, nstat)
            names(gene.weights12) <- names(statistics)
            gene.weights12[index] <- gene.weights
            index <- index2 <- rep_len(FALSE, nstat)
            names(index) <- names(index2) <- names(statistics)
            index[gene.weights12 > 0] <- TRUE
            index2[gene.weights12 < 0] <- TRUE
            gene.weights1 <- gene.weights12[gene.weights12 > 
                0]
            gene.weights2 <- gene.weights12[gene.weights12 < 
                0]
            gene.weights <- gene.weights1
        }
    }
    TWO <- !is.null(index2)
    if (is.logical(index)) 
        idx <- index
    else {
        idx <- rep_len(FALSE, nstat)
        names(idx) <- names(statistics)
        idx[index] <- TRUE
    }
    set1 <- data.frame(idx = idx, weight = NA, wt = NA)
    if (TWO) {
        if (is.logical(index2)) 
            idx2 <- index2
        else {
            idx2 <- rep_len(FALSE, nstat)
            names(idx2) <- names(statistics)
            idx2[index2] <- TRUE
        }
        set2 <- data.frame(idx = idx2, weight = NA, wt = NA)
    }
    if (length(gene.weights)) {
        set1$weight <- 0
        set1$weight[index] <- gene.weights
        set1$wt <- abs(set1$weight)/sum(abs(set1$weight))
        if (TWO) {
            set2$weight <- 0
            set2$weight[index2] <- gene.weights2
            SUM <- sum(abs(set1$weight), abs(set2$weight))
            set1$wt <- abs(set1$weight)/SUM
            set2$wt <- abs(set2$weight)/SUM
        }
    }
    ostat <- order(statistics, na.last = TRUE, decreasing = FALSE)
    statistics <- statistics[ostat]
    set1 <- set1[ostat, ]
    if (TWO) 
        set2 <- set2[ostat, ]
    n <- sum(!is.na(statistics))
    if (n == 0L) {
        message("No valid statistics")
        return(invisible())
    }
    if (n < nstat) {
        statistics <- statistics[1:n]
        set1 <- set1[1:n, ]
        if (TWO) 
            set2 <- set2[1:n, ]
    }
    r <- which(set1$idx)
    if (TWO) {
        r2 <- which(set2$idx)
        if (!length(r2)) 
            TWO <- FALSE
    }
    if (!length(r)) 
        if (TWO) {
            r <- r2
            set1 <- set2
            TWO <- FALSE
        }
        else {
            message("No selected genes to plot")
            return(invisible())
        }
    WTS <- FALSE
    wt1 <- set1$wt[r]
    len.up <- 1
    if (!anyNA(wt1)) {
        len.up <- set1$weight[r]/max(abs(set1$weight[r]))
        anydifferent <- function(x) {
            if (length(x) < 2) 
                return(FALSE)
            r <- range(x)
            (r[2] > r[1])
        }
        if (!TWO) 
            if (anydifferent(wt1)) 
                WTS <- TRUE
        if (TWO) {
            wt12 <- c(set1$wt[r], abs(set2$wt[r2]))
            if (anydifferent(wt12)) 
                WTS <- TRUE
            max.wt <- max(set1$wt[r], set2$wt[r2])
            len.up <- set1$wt[r]/max.wt
            len.down <- set2$wt[r2]/max.wt
        }
    }
    pos.dir <- all(len.up > 0)
    if (WTS) 
        shift <- 0.1
    else shift <- 0
    quantiles <- sort(quantiles)
    ALP <- alpha
    ALP <- min(ALP, 1)
    ALP <- max(ALP, 0.1)
    if (is.null(col.bars)) {
        if (TWO) {
            col.bars <- c("red", "blue")
            if (WTS) 
                col.bars.alpha <- c(rgb(1, 0, 0, alpha = ALP), 
                  rgb(0, 0, 1, alpha = ALP))
            else col.bars.alpha <- col.bars
        }
        else {
            col.bars <- "black"
            if (WTS) 
                col.bars.alpha <- rgb(0, 0, 0, alpha = ALP)
            else col.bars.alpha <- col.bars
        }
    }
    else {
        if (TWO) {
            if (length(col.bars) == 1) 
                col.bars <- rep(col.bars, 2)
            RGB <- col2rgb(col.bars)/255
            red <- RGB[1, 1]
            green <- RGB[2, 1]
            blue <- RGB[3, 1]
            red2 <- RGB[1, 2]
            green2 <- RGB[2, 2]
            blue2 <- RGB[3, 2]
            if (WTS) 
                col.bars.alpha <- c(rgb(red, green, blue, alpha = ALP), 
                  rgb(red2, green2, blue2, alpha = ALP))
            else col.bars.alpha <- col.bars
        }
        else {
            RGB <- col2rgb(col.bars)/255
            red <- RGB[1, 1]
            green <- RGB[2, 1]
            blue <- RGB[3, 1]
            if (WTS) 
                col.bars.alpha <- rgb(red, green, blue, alpha = ALP)
            else col.bars.alpha <- col.bars
        }
    }
    ylim.worm <- ylim <- c(-1, 1)
    ylab.worm <- ""
    xlab.worm <- xlab
    if (!TWO) 
        ylim.worm <- c(0, 1)
    if (worm) {
        ylim.worm <- c(-2.1, 2.1)
        if (!TWO) 
            ylim.worm <- c(0, 2.1)
    }
    ylim[2] <- ylim[2] + 0.5
    if (TWO) 
        ylim[1] <- ylim[1] - 0.5
    if (TWO) 
        plot(1:n, xlim = c(0, n), ylim = c(ylim.worm[1] - shift, 
            ylim.worm[2] + shift), type = "n", axes = FALSE, 
            xlab = xlab.worm, ylab = ylab.worm, ...)
    if (!TWO) 
        plot(1:n, xlim = c(0, n), ylim = c(ylim.worm[1] - shift * 
            (!pos.dir), ylim.worm[2] + shift * pos.dir), type = "n", 
            axes = FALSE, xlab = xlab.worm, ylab = ylab.worm, 
            ...)
    npos <- sum(statistics > quantiles[2])
    nneg <- sum(statistics < quantiles[1])
    lwd <- 50/length(r)
    lwd <- min(1.9, lwd)
    lwd <- max(0.2, lwd)
    if (TWO) {
        lwd2 <- 50/length(r2)
        lwd2 <- min(1.9, lwd2)
        lwd2 <- max(0.2, lwd2)
        lwd <- lwd2 <- min(lwd, lwd2)
    }
    barlim <- ylim[2] - c(1.5, 0.5)
    if (!pos.dir) {
        rect.yb <- 0.5
        rect.yt <- 1
        rect(nneg + 0.5, rect.yb, n - npos + 0.5, rect.yt, col = "lightgray", 
            border = NA)
        if (nneg) 
            rect(0.5, rect.yb, nneg + 0.5, rect.yt, col = rectangle1col, 
                border = NA)
        if (npos) 
            rect(n - npos + 0.5, rect.yb, n + 0.5, rect.yt, col = rectangle2col, 
                border = NA)
        segments(r, barlim[2]/2, r, barlim[2], lwd = lwd, col = col.bars.alpha[1])
        segments(r, barlim[2]/2 - shift, r, barlim[2]/2 * (1 + 
            len.up) - shift, lwd = lwd, col = col.bars[1])
    }
    if (pos.dir) {
        rect.yb <- -0.5
        if (!TWO) 
            rect.yb <- 0
        rect.yt <- 0.5
        rect(nneg + 0.5, rect.yb, n - npos + 0.5, rect.yt, col = "lightgray", 
            border = NA)
        if (nneg) 
            rect(0.5, rect.yb, nneg + 0.5, rect.yt, col = rectangle1col, 
                border = NA)
        if (npos) 
            rect(n - npos + 0.5, rect.yb, n + 0.5, rect.yt, col = rectangle2col, 
                border = NA)
        segments(r, barlim[1], r, barlim[2]/2, lwd = lwd, col = col.bars.alpha[1])
        segments(r, barlim[2]/2 + shift, r, barlim[2]/2 * (1 + 
            len.up) + shift, lwd = lwd, col = col.bars[1])
    }
    if (TWO) {
        barlim2 <- ylim[1] + c(0.5, 1.5)
        segments(r2, barlim2[1]/2, r2, barlim2[2], lwd = lwd2, 
            col = col.bars.alpha[2])
        segments(r2, barlim2[1]/2 * (1 + len.down) - shift, r2, 
            barlim2[1]/2 - shift, lwd = lwd2, col = col.bars[2])
    }
    lab.at <- 0
    if (!TWO) 
        lab.at <- 0.5
    axis(side = 2, at = lab.at, padj = 3.8, cex.axis = 0.85, 
        labels = labels[1], tick = FALSE)
    axis(side = 4, at = lab.at, padj = -3.8, cex.axis = 0.85, 
        labels = labels[2], tick = FALSE)
    prob <- (0:10)/10
    axis(at = seq(1, n, len = 11), side = 1, cex.axis = 0.7, 
        las = 2, labels = format(quantile(statistics, p = prob), 
            digits = 1))
    if (worm) {
        rescale <- function(x, newrange, oldrange = range(x)) {
            newrange[1] + (x - oldrange[1])/(oldrange[2] - oldrange[1]) * 
                (newrange[2] - newrange[1])
        }
        if (!WTS) {
            ave.enrich1 <- length(r)/n
            worm1 <- tricubeMovingAverage(set1$idx, span = span.worm)/ave.enrich1
            if (TWO) {
                ave.enrich2 <- length(r2)/n
                worm2 <- tricubeMovingAverage(-set2$idx, span = span.worm)/ave.enrich2
            }
        }
        if (WTS) {
            ave.enrich1 <- mean(set1$wt)
            worm1 <- tricubeMovingAverage(set1$wt, span = span.worm)/ave.enrich1
            if (TWO) {
                ave.enrich2 <- mean(set2$wt)
                worm2 <- tricubeMovingAverage(-set2$wt, span = span.worm)/ave.enrich2
            }
        }

        if (is.null(worm.max))
            worm.max <- max(worm1)
        r.worm1 <- c(0, worm.max)

        # max.worm1 <- max(worm1)
        # r.worm1 <- c(0, max.worm1)

        worm1.scale <- rescale(worm1, newrange = c(1.1 + shift * 
            pos.dir, 2.1 + shift * pos.dir), oldrange = r.worm1)
        if (TWO) {
            min.worm2 <- min(worm2)
            r.worm2 <- c(min.worm2, 0)
            worm2.scale <- rescale(worm2, newrange = c(-2.1 - 
                shift, -1.1 - shift), oldrange = r.worm2)
        }
        if (!TWO) {
            lines(x = 1:n, y = worm1.scale, col = col.bars[1], 
                lwd = 2)
            abline(h = rescale(1, newrange = c(1.1 + shift * 
                pos.dir, 2.1 + shift * pos.dir), oldrange = r.worm1), 
                lty = 2)
            axis(side = 2, at = c(1.1 + shift * pos.dir, 2.1 + 
                shift * pos.dir), cex.axis = 0.8, labels = c(0, 
                format(worm.max, digits = 2)))
            axis(side = 2, labels = "Enrichment", at = 1.6 + 
                shift * pos.dir, padj = -0.6, tick = FALSE, cex.axis = 0.8)
        }
        if (TWO) {
            lines(x = 1:n, y = worm1.scale, col = col.bars[1], 
                lwd = 2)
            abline(h = rescale(1, newrange = c(1.1 + shift, 2.1 + 
                shift), oldrange = r.worm1), lty = 2)
            lines(x = 1:n, y = worm2.scale, col = col.bars[2], 
                lwd = 2)
            abline(h = rescale(-1, newrange = c(-2.1 - shift, 
                -1.1 - shift), oldrange = r.worm2), lty = 2)
            axis(side = 2, at = c(1.1 + shift, 2.1 + shift), 
                cex.axis = 0.7, labels = c(0, format(worm.max, 
                  digits = 2)))
            axis(side = 2, at = c(-1.1 - shift, -2.1 - shift), 
                cex.axis = 0.7, labels = c(0, format(-min.worm2, 
                  digits = 2)))
            axis(side = 2, labels = "Enrichment", at = 1.6 + 
                shift, tick = FALSE, padj = -0.6, cex.axis = 0.7)
            axis(side = 2, labels = "Enrichment", at = -1.6 - 
                shift, tick = FALSE, padj = -0.6, cex.axis = 0.7)
        }
    }
    if (WTS) {
        if (!TWO) {
            if (pos.dir) {
                axis(side = 2, at = c(0.5 + shift, 1 + shift), 
                  cex.axis = 0.48, padj = 1.6, labels = c(0, 
                    format(max(set1$weight[r]), digits = 2)))
                axis(side = 2, labels = weights.label[1], at = 0.75 + 
                  shift, padj = 1, tick = FALSE, cex.axis = 0.5)
            }
            if (!pos.dir) {
                axis(side = 2, at = c(0 - shift, 0.5 - shift), 
                  cex.axis = 0.48, padj = 1.6, labels = c(format(min(set1$weight[r]), 
                    digits = 2), 0))
                axis(side = 2, labels = weights.label[1], at = 0.25 - 
                  shift, padj = 1, tick = FALSE, cex.axis = 0.5)
            }
        }
        if (TWO) {
            max.weight <- max(set1$weight[r], abs(set2$weight[r2]))
            axis(side = 2, at = c(0.5 + shift, 1 + shift), cex.axis = 0.43, 
                labels = c(0, format(max.weight, digits = 2, 
                  scientific = FALSE)), padj = 1.6)
            axis(side = 2, labels = weights.label[1], at = 0.75 + 
                shift, padj = 1, tick = FALSE, cex.axis = 0.46)
            axis(side = 2, at = c(-0.5 - shift, -1 - shift), 
                cex.axis = 0.43, labels = c(0, format(-max.weight, 
                  digits = 2, scientific = FALSE)), padj = 1.6)
            axis(side = 2, labels = weights.label[1], at = -0.75 - 
                shift, padj = 1, tick = FALSE, cex.axis = 0.46)
        }
    }
    invisible()
}
