############################################################
# Xenium + Morphology integration, clustering & state-wise tests
# - Builds Morphology assay from per-FOV CSVs
# - Merges 13 FOVs, attaches Morphology to your integrated object
# - UMAP on morphology metrics
# - Morphology “markers” per predicted cell state (stroma & myeloid)
# - Heatmaps with asterisk overlay for significant feature–state pairs
############################################################

## ---- Libraries ----
library(Seurat)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ArchR)   # paletteDiscrete / paletteContinuous
library(ggsci)
options(future.globals.maxSize = 90000 * 1024^2) # ~90 GB

## ---- Helpers ----
pal_disc <- function(vals, set = "stallion") paletteDiscrete(unique(vals), set = set)

safe_pdf <- function(file, width = 7, height = 5, code) {
  pdf(file, useDingbats = FALSE, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(code)
}

add_morph_assay <- function(path, disease, fov = "fov") {
  # infer CSV name: _<folder>_cell_morphology_metrics.csv (sits next to <folder>)
  folder <- basename(path)
  morph_csv <- file.path(dirname(path), paste0("_", folder, "_cell_morphology_metrics.csv"))
  if (!file.exists(morph_csv)) stop("Missing morphology CSV: ", morph_csv)

  obj <- LoadXenium(path, fov = fov)
  # morphology file: rows = metrics, cols = cells (after transpose)
  morp <- t(read.csv(morph_csv, row.names = 1, check.names = FALSE))
  obj[["Morphology"]] <- CreateAssayObject(counts = as.matrix(morp))
  obj@assays$Morphology@key <- "morph_"
  obj <- subset(obj, subset = nCount_Xenium > 0)
  obj$disease <- disease
  obj
}

## ---- Load all FOVs with morphology ----
objs <- list(
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_1__20240621__170414", "Mild"),
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_3__20240621__170414", "Moderate"),
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_5__20240621__170414", "Mild"),
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_6__20240621__170414", "Mild"),

  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_1__20240621__170414", "Moderate"),
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_2__20240621__170414", "Mild"),
  add_morph_assay("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_3__20240621__170414", "Mild"),

  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_1__20240628__212431", "Severe"),
  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_2__20240628__212431", "Severe"),
  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_3__20240628__212431", "Severe"),

  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_1__20240628__212431", "Moderate"),
  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_3__20240628__212431", "Moderate"),
  add_morph_assay("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_5__20240628__212431", "Mild")
)

xenium.obj.m <- Reduce(function(a,b) merge(a, y = b), objs)
xenium.obj.m <- subset(xenium.obj.m, subset = nCount_Xenium > 10)

## ---- Attach Morphology to your integrated object ----
xenium.obj <- readRDS("../analysis_amrute/xenium.obj.integrated.rds")
xenium.obj[["Morphology"]] <- xenium.obj.m[["Morphology"]]
saveRDS(xenium.obj, "xenium.obj.integrated.morphology.rds")

## ---- Morphology UMAP (global) ----
DefaultAssay(xenium.obj) <- "Morphology"
xenium.obj <- ScaleData(xenium.obj, features = rownames(xenium.obj[["Morphology"]]), verbose = FALSE)
xenium.obj <- RunPCA(
  xenium.obj,
  assay = "Morphology",
  features = rownames(xenium.obj[["Morphology"]]),
  reduction.name = "pca_morphology", verbose = FALSE
)
xenium.obj <- RunUMAP(xenium.obj, reduction = "pca_morphology", dims = 1:10, assay = "Morphology")

DimPlot(xenium.obj, group.by = "predicted.id", reduction = "umap") +
  theme(aspect.ratio = 1)

## ---- Load compartment-specific mapped objects (for labels) ----
xenium.Myeloid <- readRDS("../analysis_amrute/xenium.Myeloid.rds")
xenium.stroma  <- readRDS("../analysis_amrute/xenium.stroma.rds")

############################################################
# STROMA — Morphology markers per mapped cell.state
############################################################
Idents(xenium.obj) <- "predicted.id"
stroma <- subset(xenium.obj, idents = c("Fibroblast2","Fibroblast1","SMCPericyte","ModSMC"))

# bring mapped cell.state labels (barcode intersection assumed)
stroma$cell.state <- xenium.stroma$predicted.celltype[Cells(stroma)]
Idents(stroma) <- "cell.state"

# Find morphology “markers” on scaled data
morph_markers_stroma <- FindAllMarkers(
  stroma,
  assay = "Morphology",
  slot  = "scale.data",
  only.pos = FALSE
)

morph_markers_stroma <- subset(morph_markers_stroma, p_val_adj < 0.05)

# Average morphology per cell.state
DefaultAssay(stroma) <- "Morphology"
avg_stroma <- AverageExpression(
  stroma,
  assays   = "Morphology",
  features = rownames(stroma[["Morphology"]]),
  group.by = "cell.state",
  slot     = "counts"   # use raw morphology values for display
)[[1]]

# Make an asterisk overlay for significant (feature, state)
morph_markers_stroma$gene <- sub("\\.\\d+$", "", morph_markers_stroma$gene)
rownames(avg_stroma)     <- sub("\\.\\d+$", "", rownames(avg_stroma))

anno_stroma <- matrix(
  "", nrow = nrow(avg_stroma), ncol = ncol(avg_stroma),
  dimnames = list(rownames(avg_stroma), colnames(avg_stroma))
)
if (nrow(morph_markers_stroma) > 0) {
  for (i in seq_len(nrow(morph_markers_stroma))) {
    g <- morph_markers_stroma$gene[i]
    cl <- as.character(morph_markers_stroma$cluster[i])
    if (g %in% rownames(anno_stroma) && cl %in% colnames(anno_stroma)) {
      anno_stroma[g, cl] <- "*"
    }
  }
}

safe_pdf("./stroma_cell_states_morphology_metrics.pdf", 3.5, 2.5, {
  pheatmap::pheatmap(
    avg_stroma,
    scale          = "row",
    col            = paletteContinuous("solarExtra"),
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    fontsize_row   = 6,
    fontsize_col   = 6,
    border_color   = NA,
    legend         = TRUE,
    display_numbers= anno_stroma,
    number_color   = "black"
  )
})

############################################################
# MYELOID — Morphology markers per mapped cell.state
############################################################
myeloid <- subset(xenium.obj, idents = "Myeloid")
myeloid$cell.state <- xenium.Myeloid$predicted.celltype[Cells(myeloid)]
Idents(myeloid) <- "cell.state"

morph_markers_my <- FindAllMarkers(
  myeloid,
  assay = "Morphology",
  slot  = "scale.data",
  only.pos = FALSE
)
morph_markers_my <- subset(morph_markers_my, p_val_adj < 0.05)

# Dimension reduction on morphology (myeloid only)
myeloid <- RunPCA(myeloid, assay = "Morphology",
                  features = rownames(myeloid[["Morphology"]]),
                  reduction.name = "pca_morphology", verbose = FALSE)
myeloid <- RunUMAP(myeloid, reduction = "pca_morphology", dims = 1:10, assay = "Morphology")

# Average & asterisk overlay
DefaultAssay(myeloid) <- "Morphology"
avg_my <- AverageExpression(
  myeloid,
  assays   = "Morphology",
  features = rownames(myeloid[["Morphology"]]),
  group.by = "cell.state",
  slot     = "counts"
)[[1]]

morph_markers_my$gene <- sub("\\.\\d+$", "", morph_markers_my$gene)
rownames(avg_my)      <- sub("\\.\\d+$", "", rownames(avg_my))

anno_my <- matrix(
  "", nrow = nrow(avg_my), ncol = ncol(avg_my),
  dimnames = list(rownames(avg_my), colnames(avg_my))
)
if (nrow(morph_markers_my) > 0) {
  for (i in seq_len(nrow(morph_markers_my))) {
    g <- morph_markers_my$gene[i]
    cl <- as.character(morph_markers_my$cluster[i])
    if (g %in% rownames(anno_my) && cl %in% colnames(anno_my)) {
      anno_my[g, cl] <- "*"
    }
  }
}

safe_pdf("./myeloid_cell_states_morphology_metrics.pdf", 3.5, 2.5, {
  pheatmap::pheatmap(
    avg_my,
    scale          = "row",
    col            = paletteContinuous("solarExtra"),
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    fontsize_row   = 6,
    fontsize_col   = 6,
    border_color   = NA,
    legend         = TRUE,
    display_numbers= anno_my,
    number_color   = "black"
  )
})
