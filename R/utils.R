################################################################################
# Shared helper functions used across the analysis scripts
#
# Paper : Amrute et al., Science (2026) | doi:10.1126/science.adx1736
#
# Usage : source(file.path(repo_dir, "R", "utils.R"))
#         where `repo_dir` points to the root of this repository.
#
# These functions replace blocks that were previously copy-pasted in many
# scripts. Each one reproduces the original computation exactly.
################################################################################

suppressPackageStartupMessages(library(ggplot2))

#' Gene-set z-score
#'
#' Scales each gene across cells (z-score), sets genes with zero variance to 0,
#' and averages the z-scores of the matched genes for every cell. Gene names are
#' matched case-insensitively, so one list works for human and mouse symbols.
#'
#' @param expdata Gene x cell expression matrix, e.g. GetAssayData(obj) for the
#'   current default assay.
#' @param genes Character vector of gene symbols.
#' @return Named numeric vector (one score per cell).
#' @examples
#' DefaultAssay(obj) <- "RNA"
#' obj$FOAM <- gene_set_zscore(GetAssayData(obj), c("FABP4", "FABP5", "GPNMB"))
gene_set_zscore <- function(expdata, genes) {
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  if (length(zz) == 0) stop("None of the genes were found in the expression matrix.")
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  colSums(geneExp) / length(zz)
}

#' Map cluster IDs to annotation labels
#'
#' @param clusters Vector of cluster IDs (factor, character or numeric), e.g.
#'   obj$SCT_snn_res.0.2. Names (cell barcodes) are preserved.
#' @param labels Named character vector: names are cluster IDs, values are labels,
#'   e.g. c("0" = "TCells", "1" = "Myeloid").
#' @return Character vector of labels, same length and names as `clusters`.
annotate_clusters <- function(clusters, labels) {
  ids <- as.character(clusters)
  missing <- setdiff(unique(ids), names(labels))
  if (length(missing) > 0) {
    stop("No label for cluster(s): ", paste(missing, collapse = ", "))
  }
  stats::setNames(unname(labels[ids]), names(clusters))
}

#' Stacked proportion bar plot (cell composition)
#'
#' @param meta Data frame of cell metadata, e.g. obj@meta.data.
#' @param x Column used for the x axis (e.g. "condition", "sampleID").
#' @param fill Column used for the fill (e.g. "cell.type").
#' @param colors Fill colours passed to scale_fill_manual().
#' @return A ggplot object.
plot_composition <- function(meta, x, fill, colors) {
  ggplot(meta, aes(x = !!rlang::sym(x), fill = !!rlang::sym(fill))) +
    geom_bar(position = "fill") + theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_manual(values = colors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

#' Read a condition DE table (FindAllMarkers output) and flag significance
#'
#' @param path CSV written by FindAllMarkers().
#' @param cell Cell-type label added as column `cell`.
#' @param lfc_cutoff Minimum avg_log2FC for "Significant" (with p_val_adj < 0.05).
#' @return Data frame with added columns `cell`, `sigpvalue`, `sig`.
read_condition_de <- function(path, cell, lfc_cutoff = 0.25) {
  de <- readr::read_delim(path, ",", show_col_types = FALSE)
  de$cell <- cell
  de$sigpvalue <- ifelse(de$p_val_adj < 0.05, "p < 0.05", "p > 0.05")
  de$sig <- ifelse(de$p_val_adj < 0.05 & de$avg_log2FC > lfc_cutoff,
                   "Significant", "Not Significant")
  de
}
