################################################################################
# Mouse-to-human mapping: Wirka et al. 2019 SMC lineage-traced atherosclerosis scRNA-seq
#
# Paper : Amrute et al., Science (2026) | doi:10.1126/science.adx1736
# Part  : Origin of FAP+ modulated SMCs (mouse lineage tracing)
#
# Purpose
#   Re-processes GSE131776 (Myh11-lineage-traced SMCs; wild-type and SMC-specific
#   Tcf21 knockout; baseline, 8 and 16 weeks of high-fat diet), converts mouse
#   genes to human orthologs, and projects SMC-derived cells onto the human
#   stromal CITE-seq reference to follow which human states they occupy.
#
# Inputs
#   GSE131776_mouse_scRNAseq_wirka_et_al_GEO.txt  (GEO GSE131776)
#   smc_fib_annotated.rds  (human stromal reference; needs `spca` + `rna.umap` model)
#
# Outputs
#   TQ_Wirka_mapped.rds
#
# Run order
#   Upstream  : human_citeseq/06_stroma_cell_states_FAP.R
#   Downstream: none
#
# Notes
#   Interactive analysis script: run section by section (e.g. in RStudio).
################################################################################

## ---- Shared helpers (R/utils.R) ----
repo_dir <- "."  # EDIT: path to the root of this repository
source(file.path(repo_dir, "R", "utils.R"))

## ---- Libraries (unique + sufficient) ----
library(Seurat)        # core single-cell workflow
library(ggplot2)       # plotting
library(patchwork)     # plot composition
library(ggpubr)        # publication-friendly plots
library(dplyr)         # data wrangling
library(sctransform)   # not used directly
library(pheatmap)      # heatmaps
library(Matrix)        # sparse matrices
library(RColorBrewer)  # color palettes
library(scales)        # scaling helpers (limits/transform)
library(data.table)    # fast table ops (optional)
library(stats)         # base stats
library(Nebulosa)      # density plots (plot_density)
library(ggsci)         # extra palettes (optional)
library(ArchR)         # palettes: paletteDiscrete/paletteContinuous
library(biomaRt)       # not used directly
library(nichenetr)     # convert_mouse_to_human_symbols()

## ==========================================================
## 1) Read Wirka raw counts and create Seurat object
## ==========================================================
# Expect a gene-by-cell matrix in TSV format
raw_counts <- read.table(file = "./GSE131776_mouse_scRNAseq_wirka_et_al_GEO.txt",
                         sep = "\t")
head(raw_counts)

# Create object and basic QC thresholds
mydata <- CreateSeuratObject(counts = raw_counts, min.cells = 5, min.features = 500)

# Compute mitochondrial % (mouse genes typically "^mt-")
mydata <- PercentageFeatureSet(mydata, pattern = "^mt-", col.name = "percent.mt")

# Filter by mt%, feature count
mydata <- subset(mydata, percent.mt < 7.5 & nFeature_RNA > 500 & nFeature_RNA < 3500)

# Extract sample code from column names (suffix after last dot)
mydata$sample <- sub(".*\\.", "", colnames(mydata))
unique(mydata$sample)

## ==========================================================
## 2) Map sample codes to experimental conditions; subset to SMC samples
## ==========================================================
# Map sample integers ("1".."18") to conditions
mydata_condition_labels <- c(
  "1" = "wt_SMC_baseline", "2" = "wt_nonSMC_baseline", "3" = "wt_SMC_baseline", "4" = "wt_nonSMC_baseline",
  "5" = "wt_SMC_8wk", "6" = "wt_nonSMC_8wk", "7" = "wt_SMC_8wk", "8" = "wt_nonSMC_8wk",
  "9" = "ko_SMC_8wk", "10" = "ko_nonSMC_8wk", "11" = "wt_SMC_16wk", "12" = "wt_nonSMC_16wk",
  "13" = "ko_SMC_16wk", "14" = "ko_nonSMC_16wk", "15" = "wt_SMC_16wk", "16" = "ko_SMC_16wk",
  "17" = "ko_SMC_16wk", "18" = "ko_nonSMC_16wk"
)
mydata$condition <- annotate_clusters(mydata$sample, mydata_condition_labels)

# Keep SMC lineage-traced samples only
Idents(mydata) <- "condition"
mydata <- subset(mydata, idents = c("wt_SMC_baseline","wt_SMC_8wk","ko_SMC_8wk","wt_SMC_16wk","ko_SMC_16wk"))

## ==========================================================
## 3) Load human reference (SMC/Fib annotated) and update object
## ==========================================================
# Reference is expected to contain RNA data + SPCA + RNA UMAP
sample <- readRDS("./smc_fib_annotated.rds")

# Update mydata to latest Seurat structure
mydata <- UpdateSeuratObject(mydata)

## ==========================================================
## 4) Convert mouse gene symbols to human, rebuild object with same metadata
## ==========================================================
mouse_rna_matrix <- mydata@assays[["RNA"]]@counts
# convert_mouse_to_human_symbols() is provided by nichenetr
rownames(mouse_rna_matrix) <- mouse_rna_matrix %>%
  rownames() %>%
  convert_mouse_to_human_symbols()

# Drop rows/cols that became NA during conversion
mouse_rna_matrix <- mouse_rna_matrix %>% .[!is.na(rownames(mouse_rna_matrix)),
                                            !is.na(colnames(mouse_rna_matrix))]

# Reconstruct Seurat object and carry over metadata
mouse_coronary_new <- CreateSeuratObject(counts = mouse_rna_matrix)
mouse_coronary_new <- AddMetaData(object = mouse_coronary_new, metadata = mydata@meta.data)

# Intersect features with reference to ensure consistent space
DefaultAssay(mouse_coronary_new) <- "RNA"
DefaultAssay(sample) <- "RNA"
mouse_coronary_new <- subset(mouse_coronary_new, features = rownames(sample))

## ==========================================================
## 5) Light preprocessing before mapping
## ==========================================================
DefaultAssay(mouse_coronary_new) <- "RNA"
mouse_coronary_new <- NormalizeData(mouse_coronary_new)
mouse_coronary_new <- FindVariableFeatures(mouse_coronary_new, selection.method = "vst", nfeatures = 3000)
mouse_coronary_new <- ScaleData(mouse_coronary_new)

## ==========================================================
## 6) Reference mapping to human (Seurat v5 style)
## ==========================================================
human_coronary <- sample

anchors <- FindTransferAnchors(
  reference = human_coronary,
  query = mouse_coronary_new,
  normalization.method = "LogNormalize",
  reference.reduction = "spca",
  dims = 1:50,
  reference.assay = "RNA"
)

mouse_coronary_new <- MapQuery(
  anchorset = anchors,
  query = mouse_coronary_new,
  reference = human_coronary,
  refdata = list(
    celltype = "cell.state",   # metadata column in reference
    predicted_ADT = "ADT"      # optional flow-through if present
  ),
  reference.reduction = "spca",
  reduction.model = "rna.umap" # must exist in the reference
)

## ==========================================================
## 7) Factor ordering, subset to WT groups for some plots
## ==========================================================
mouse_coronary_new$condition <- factor(
  mouse_coronary_new$condition,
  levels = c("wt_SMC_baseline","wt_SMC_8wk","wt_SMC_16wk","ko_SMC_8wk","ko_SMC_16wk")
)

Idents(mouse_coronary_new) <- "condition"
mouse_coronary_new <- subset(mouse_coronary_new, idents = c("wt_SMC_baseline","wt_SMC_8wk","wt_SMC_16wk"))

## ==========================================================
## 8) Plots: predicted cell types, composition, features
## ==========================================================
# UMAP labeled by predicted reference cell types
DimPlot(
  mouse_coronary_new, reduction = "ref.umap", group.by = "predicted.celltype",
  label.size = 4, cols = paletteDiscrete(unique(sample$cell.state), set = "stallion"),
  label = FALSE
)

# Composition per condition
plot_composition(mouse_coronary_new@meta.data, x = "condition", fill = "predicted.celltype",
                 colors = paletteDiscrete(unique(sample$cell.state), set = "stallion"))

# Example feature summaries
DotPlot(mouse_coronary_new, features = "FAP", group.by = "condition") + RotatedAxis()
RidgePlot(mouse_coronary_new, features = "predicted.celltype.score")

# Save mapped object
saveRDS(mouse_coronary_new, "TQ_Wirka_mapped.rds")

