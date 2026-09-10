############################################################
# Mouse (Owens: PMID: 32674599) → Human (CITE-seq) mapping
# - Converts mouse symbols → human
# - Normalizes & prepares mouse data
# - Ensures human reference has SPCA + UMAP(model) for MapQuery
# - Transfers labels, visualizes, and saves outputs
############################################################

## Packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)
library(pheatmap)
library(Matrix)
library(RColorBrewer)
library(scales)
library(data.table)
library(stats)
library(ggsci)
library(ArchR)     # for paletteDiscrete / paletteContinuous
library(biomaRt)
library(nichenetr)
library(tidyverse)

## ---------- IO ----------
human_ref_rds <- "./smc_fib_annotated.rds"
mouse_rds     <- "./healthy.disease.OwensLab.only.making.final.graphs15dim.new.idents.v2.rds"
out_rds       <- "Owens_ApoE_plaque_mapped_toCITEseqHuman.rds"

## ---------- Load data ----------
sample <- readRDS(human_ref_rds)           # human CITE-seq stroma reference
mydata <- readRDS(mouse_rds) %>% UpdateSeuratObject()

# Keep the Owens origins of interest
Idents(mydata) <- "origin"
mydata <- subset(mydata, idents = c("SMC_PDGF_WT_eYFP_Positive","tdTomato","SMC_Klf4_WT_eYFP_Positive","GFP"))

## ---------- Mouse → Human gene symbol conversion ----------
# Prefer nichenetr's converter; fall back to biomaRt if needed
convert_m2h <- function(genes_mouse) {
  if (exists("convert_mouse_to_human_symbols")) {
    out <- convert_mouse_to_human_symbols(genes_mouse)
  } else {
    message("Using biomaRt fallback for mouse→human symbol conversion")
    m <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    h <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    map <- getLDS(attributes = c("mgi_symbol"),
                  filters    = "mgi_symbol",
                  values     = genes_mouse,
                  mart       = m,
                  attributesL= c("hgnc_symbol"),
                  martL      = h,
                  uniqueRows = TRUE)
    colnames(map) <- c("mouse","human")
    out <- setNames(map$human, map$mouse)[genes_mouse]
  }
  # Ensure unique, upper-case standard human symbols where possible
  out <- toupper(out)
  return(out)
}

mouse_rna_matrix <- mydata[["RNA"]]@counts
rownames(mouse_rna_matrix) <- convert_m2h(rownames(mouse_rna_matrix))

# Drop NAs created by unmapped symbols
mouse_rna_matrix <- mouse_rna_matrix[!is.na(rownames(mouse_rna_matrix)), , drop = FALSE]
mouse_rna_matrix <- mouse_rna_matrix[rownames(mouse_rna_matrix) != "", , drop = FALSE]

# Build Seurat object from converted matrix and carry metadata
mouse_coronary_new <- CreateSeuratObject(counts = mouse_rna_matrix, meta.data = mydata@meta.data)

## ---------- Align gene space to human reference ----------
DefaultAssay(mouse_coronary_new) <- "RNA"
DefaultAssay(sample)             <- "RNA"
common_genes <- intersect(rownames(sample), rownames(mouse_coronary_new))
mouse_coronary_new <- subset(mouse_coronary_new, features = common_genes)

## ---------- Process the mouse dataset ----------
mouse_coronary_new <- NormalizeData(mouse_coronary_new)
mouse_coronary_new <- FindVariableFeatures(mouse_coronary_new, selection.method = "vst", nfeatures = 3000)
mouse_coronary_new <- ScaleData(mouse_coronary_new, features = rownames(mouse_coronary_new))

## ---------- Ensure human reference has SPCA + UMAP(model) ----------
# SPCA on human ref if missing
if (!"spca" %in% names(sample@reductions)) {
  message("Running RunSPCA on human reference...")
  sample <- RunSPCA(sample, assay = "RNA", verbose = FALSE)
}

# Ensure a UMAP model exists for MapQuery (reduction.model = "rna.umap")
if (!"rna.umap" %in% names(sample@reductions) ||
    is.null(sample@reductions[["rna.umap"]]@misc$model)) {
  message("Running UMAP with return.model=TRUE on the reference...")
  red_use <- if ("spca" %in% names(sample@reductions)) "spca" else "pca"
  sample <- RunUMAP(sample, reduction = red_use, dims = 1:50,
                    reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
                    return.model = TRUE, verbose = FALSE)
}

human_coronary <- sample  # alias for clarity

## ---------- Transfer anchors & mapping ----------
anchors <- FindTransferAnchors(
  reference            = human_coronary,
  query                = mouse_coronary_new,
  normalization.method = "LogNormalize",
  reference.reduction  = "spca",
  dims                 = 1:50,
  reference.assay      = "RNA"
)

mouse_coronary_new <- MapQuery(
  anchorset           = anchors,
  query               = mouse_coronary_new,
  reference           = human_coronary,
  refdata             = list(
    celltype       = "cell.state",
    predicted_ADT  = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model     = "rna.umap"
)

## ---------- Plots ----------
# UMAP in reference space, colored by predicted cell type
p1 <- DimPlot(
  mouse_coronary_new,
  reduction = "ref.umap",
  group.by  = "predicted.celltype",
  label     = FALSE
) + scale_color_manual(values = paletteDiscrete(unique(human_coronary$cell.state), set = "stallion")) +
    ggtitle("Mouse cells mapped to human reference (predicted cell.state)")

# Composition across OWENS 'origin'
p2 <- ggplot(mouse_coronary_new@meta.data,
             aes(x = origin, fill = predicted.celltype)) +
  geom_bar(position = "fill") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.line   = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_fill_manual(values = paletteDiscrete(unique(human_coronary$cell.state), set = "stallion")) +
  ylab("Fraction") + xlab("Origin") +
  ggtitle("Predicted cell.state by Owens origin")

# Save plots
ggsave("owens_to_human_predicted_celltypes_refUMAP.png", p1, width = 7, height = 5, dpi = 300)
ggsave("owens_to_human_predicted_by_origin.png",       p2, width = 7, height = 5, dpi = 300)

# Ridge of prediction confidence
ridge <- RidgePlot(mouse_coronary_new, features = "predicted.celltype.score", ncol = 1) +
  ggtitle("Prediction score (Seurat label transfer)")
ggsave("owens_to_human_prediction_score_ridge.png", ridge, width = 7, height = 5, dpi = 300)

## ---------- Save ----------
saveRDS(mouse_coronary_new, out_rds)

## ---------- Session recap ----------
message("Mapped object saved to: ", out_rds)
message("Cells: ", ncol(mouse_coronary_new), " | Genes (common): ", length(common_genes))
if (!is.null(mouse_coronary_new$predicted.celltype)) {
  message("Predicted cell types (head): ", paste(head(unique(mouse_coronary_new$predicted.celltype)), collapse = ", "))
}
