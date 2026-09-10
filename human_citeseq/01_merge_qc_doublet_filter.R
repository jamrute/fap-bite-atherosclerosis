################################################################################
# Human coronary CITE-seq: merge samples, QC filtering and doublet removal
#
# Paper : Amrute et al., Science (2026) | doi:10.1126/science.adx1736
# Part  : Human coronary artery CITE-seq atlas
#
# Purpose
#   Loads Cell Ranger gene-expression + antibody-capture (ADT) matrices for
#   27 coronary artery samples, merges them, applies cell-level QC thresholds,
#   drops one low-quality sample (sample21), exports the object for Scrublet,
#   then adds Scrublet scores and removes likely doublets (score >= 0.25).
#
# Inputs
#   Cell Ranger output, one folder per sample under `base_cr` (GEO GSE314596)
#   ./scrublet/scrublet-scores/all.csv  (from 02_scrublet_doublet_scores.ipynb)
#
# Outputs
#   merged_preQC.rds
#   merged_postQC.h5Seurat / merged_postQC.h5ad  (input to Scrublet)
#   merged_postQC_doubletRemoval.rds
#
# Run order
#   Upstream  : none (first step)
#   Downstream: 02_scrublet_doublet_scores.ipynb, then normalization/clustering -> 03
#
# Notes
#   Interactive analysis script: run section by section (e.g. in RStudio).
#   Two passes: run up to the h5ad export, run notebook 02, then continue
#   from the Scrublet section.
################################################################################

## --- Libraries ---
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(SeuratDisk)
})

## --- Paths & sample lists ---
# EDIT: directory containing one Cell Ranger output folder per sample
base_cr <- "data/human_citeseq/cellranger"

# All sample IDs to load (add/remove here as needed)
sample_ids <- c(
  "sample10","sample11","sample14","sample16","sample18","sample19","sample20",
  "sample21","sample22","sample23","sample24","sample25","sample26","sample31",
  "sample35","sample36","sample37","sample38","sample40","sample43","sample44",
  "sample46","sample47","sample48","sample50","WashU1","WashU2"
)

## --- Helper to load 10x GEX + ADT into a Seurat object ---
load_citeseq <- function(sample_id, base_dir = base_cr) {
  dir_path <- file.path(base_dir, sample_id)
  if (!dir.exists(dir_path)) stop(paste("Missing directory:", dir_path))

  # Read10X returns a named list for multi-modal 10x outputs (e.g., "Gene Expression" and "Antibody Capture")
  x <- Read10X(data.dir = dir_path)

  # Create Seurat object from RNA; ADT is optional
  if (!"Gene Expression" %in% names(x)) stop(paste("No 'Gene Expression' matrix for", sample_id))
  obj <- CreateSeuratObject(counts = x[["Gene Expression"]], project = sample_id)
  if ("Antibody Capture" %in% names(x)) {
    obj[["ADT"]] <- CreateAssayObject(counts = x[["Antibody Capture"]])
  }

  obj$sampleID <- sample_id
  obj
}

## --- Load all samples (with informative messages) ---
message("Loading samples...")
objs <- lapply(sample_ids, function(sid) {
  message("  - ", sid)
  load_citeseq(sid)
})

## --- Merge all samples ---
# Use the first object as anchor, merge the rest. This is memory-friendlier than Reduce for large lists.
sample <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = sample_ids # prefixes barcodes with sample IDs (helps uniqueness & joins)
)

# Checkpoint: merged, unfiltered object
saveRDS(sample, file = "./merged_preQC.rds")

## --- Re-load if desired ---
# sample <- readRDS("./merged_preQC.rds")

## --- Basic QC metrics ---
sample[["percent.mt"]]   <- PercentageFeatureSet(sample, pattern = "^MT-")
sample[["percent.ribo"]] <- PercentageFeatureSet(sample, pattern = "^RP[SL]")

# Quick glance
print(sample)

VlnPlot(sample, features = c("nCount_RNA", "nCount_ADT", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sampleID")

## --- Filtering (parametrized) ---
qc_min_features   <- 200
qc_max_counts_rna <- 25000
qc_min_counts_rna <- 500
qc_max_counts_adt <- 10000
qc_max_mt         <- 20

sample_filtered <- subset(
  sample,
  subset = nFeature_RNA > qc_min_features &
           nCount_RNA   > qc_min_counts_rna &
           nCount_RNA   < qc_max_counts_rna &
           (is.na(nCount_ADT) | nCount_ADT < qc_max_counts_adt) &
           percent.mt   < qc_max_mt
)

VlnPlot(sample_filtered, features = c("nCount_RNA", "nCount_ADT", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "sampleID")

# Per-sample distributions (no dots, no legend)
for (feat in c("nFeature_RNA","nCount_RNA","nFeature_ADT","nCount_ADT","percent.mt","percent.ribo")) {
  if (feat %in% colnames(sample_filtered@meta.data)) {
    print(VlnPlot(sample_filtered, features = feat, group.by = "sampleID", pt.size = 0) + NoLegend())
  }
}

## --- Drop poor-quality samples (editable list) ---
drop_samples <- c("sample21")  # add others here as needed
keep_samples <- setdiff(unique(sample_filtered$sampleID), drop_samples)

Idents(sample_filtered) <- "sampleID"
sample_filtered <- subset(sample_filtered, idents = keep_samples)
print(sample_filtered)

## --- Save as h5Seurat (+ optional h5ad) ---
SaveH5Seurat(sample_filtered, filename = "./merged_postQC.h5Seurat", overwrite = TRUE)
Convert("./merged_postQC.h5Seurat", dest = "h5ad", overwrite = TRUE)

## --- Pass 2: run 02_scrublet_doublet_scores.ipynb on merged_postQC.h5ad first ---
## --- Scrublet (doublet) integration: join by cell barcode ---
# Expect a CSV with rownames or a "barcode" column matching colnames(sample_filtered)
# Read Scrublet scores and expose barcodes as an explicit column
scrub <- read.csv("./scrublet/scrublet-scores/all.csv", header = TRUE, row.names = 1, check.names = FALSE)
scrub$barcode <- rownames(scrub)

# Cell names carry sample prefixes (add.cell.ids above); Scrublet barcodes must match them.
md <- sample_filtered@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  left_join(scrub %>% select(barcode, scrublet_score, scrublet_cluster_score, bh_pval),
            by = "barcode")

# Write back with ordering preserved
stopifnot(nrow(md) == ncol(sample_filtered))
sample_filtered@meta.data <- md %>%
  select(-barcode) %>%
  as.data.frame()

# Sanity plots
if ("scrublet_score" %in% colnames(sample_filtered@meta.data)) {
  print(VlnPlot(sample_filtered, group.by = "sampleID", features = "scrublet_score", pt.size = 0) + NoLegend())
}

## --- Doublet filtering (parametrized) ---
doublet_max <- 0.25
sample_filtered2 <- subset(sample_filtered, subset = is.na(scrublet_score) | scrublet_score < doublet_max)

# Quick checks post-doublet removal
for (feat in c("nFeature_RNA","nCount_RNA","nFeature_ADT","nCount_ADT","percent.mt","percent.ribo","scrublet_score")) {
  if (feat %in% colnames(sample_filtered2@meta.data)) {
    print(VlnPlot(sample_filtered2, features = feat, group.by = "sampleID", pt.size = 0) + NoLegend())
  }
}

saveRDS(sample_filtered2, "./merged_postQC_doubletRemoval.rds")
