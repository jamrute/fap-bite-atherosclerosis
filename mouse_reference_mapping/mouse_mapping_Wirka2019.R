############################################################
# Mouse → Human mapping (Wirka et al.) and visualization
# - Deduplicated libraries
# - Original logic preserved
# - Clear, step-by-step comments
# Notes/assumptions:
#   * `convert_mouse_to_human_symbols()` exists in your env.
#   * `paletteDiscrete()` / `paletteContinuous()` exist (e.g., from ArchR utils).
#   * Reference object `sample` has SPCA + RNA UMAP (as used by MapQuery).
############################################################

## ---- Libraries (unique + sufficient) ----
library(Seurat)        # core single-cell workflow
library(ggplot2)       # plotting
library(patchwork)     # plot composition
library(ggpubr)        # publication-friendly plots
library(dplyr)         # data wrangling
library(sctransform)   # (not used directly here but kept if needed)
library(pheatmap)      # heatmaps
library(Matrix)        # sparse matrices
library(RColorBrewer)  # color palettes
library(scales)        # scaling helpers (limits/transform)
library(data.table)    # fast table ops (optional)
library(stats)         # base stats
library(Nebulosa)      # density plots (plot_density)
library(ggsci)         # extra palettes (optional)
library(ArchR)         # palettes: paletteDiscrete/paletteContinuous
library(biomaRt)       # (not directly used; kept if gene ID mapping needed)

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
fun <- function(x) {
  if (x == "1")  {"wt_SMC_baseline"}
  else if (x == "2")  {"wt_nonSMC_baseline"}
  else if (x == "3")  {"wt_SMC_baseline"}
  else if (x == "4")  {"wt_nonSMC_baseline"}
  else if (x == "5")  {"wt_SMC_8wk"}
  else if (x == "6")  {"wt_nonSMC_8wk"}
  else if (x == "7")  {"wt_SMC_8wk"}
  else if (x == "8")  {"wt_nonSMC_8wk"}
  else if (x == "9")  {"ko_SMC_8wk"}
  else if (x == "10") {"ko_nonSMC_8wk"}
  else if (x == "11") {"wt_SMC_16wk"}
  else if (x == "12") {"wt_nonSMC_16wk"}
  else if (x == "13") {"ko_SMC_16wk"}
  else if (x == "14") {"ko_nonSMC_16wk"}
  else if (x == "15") {"wt_SMC_16wk"}
  else if (x == "16") {"ko_SMC_16wk"}
  else if (x == "17") {"ko_SMC_16wk"}
  else if (x == "18") {"ko_nonSMC_16wk"}
}
mydata$condition <- mapply(fun, mydata$sample)

# Keep SMC-related groups (as in your code)
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
# Assumes helper: convert_mouse_to_human_symbols()
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
ggplot(mouse_coronary_new@meta.data,
       aes(x = condition, fill = predicted.celltype)) +
  geom_bar(position = "fill") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = paletteDiscrete(unique(sample$cell.state), set = "stallion")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Example feature summaries
DotPlot(mouse_coronary_new, features = "FAP", group.by = "condition") + RotatedAxis()
RidgePlot(mouse_coronary_new, features = "predicted.celltype.score")

# Save mapped object
saveRDS(mouse_coronary_new, "TQ_Wirka_mapped.rds")

