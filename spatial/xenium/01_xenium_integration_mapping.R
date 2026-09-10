############################################################
# Xenium multi-fov integration, clustering, mapping & visuals
# - Libraries deduplicated & ordered
# - Fixes undefined objects (xenium.obj2 → xenium.obj)
# - Guards for missing references; consistent palettes
# - Keeps your original workflow & outputs
############################################################

## ---- Libraries ----
library(Seurat)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ArchR)     # paletteDiscrete / paletteContinuous
library(ggsci)
options(future.globals.maxSize = 90000 * 1024^2) # ~90 GB

## ---- Helpers ----
pal_disc <- function(vals, set = "stallion") paletteDiscrete(unique(vals), set = set)

safe_pdf <- function(file, width = 7, height = 5, code) {
  pdf(file, useDingbats = FALSE, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(code)
}

load_xenium <- function(path, disease, fov = "fov") {
  obj <- LoadXenium(path, fov = fov)
  subset(obj, subset = nCount_Xenium > 0) |>
    AddMetaData(metadata = setNames(list(disease), "disease"))
}

## ---- Load all fovs ----
s1  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_1__20240621__170414", "Mild")
s2  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_3__20240621__170414", "Moderate")
s3  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_5__20240621__170414", "Mild")
s4  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0005092__Region_6__20240621__170414", "Mild")

s5  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_1__20240621__170414", "Moderate")
s6  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_2__20240621__170414", "Mild")
s7  <- load_xenium("../xenium_data/SR005130_Lavine1_20240621_Xenium/output-XETG00113__0010699__Region_3__20240621__170414", "Mild")

s8  <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_1__20240628__212431", "Severe")
s9  <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_2__20240628__212431", "Severe")
s10 <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0005113__Region_3__20240628__212431", "Severe")

s11 <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_1__20240628__212431", "Moderate")
s12 <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_3__20240628__212431", "Moderate")
s13 <- load_xenium("../xenium_data/SR005130_Lavine2_20240628_Xenium/output-XETG00113__0010641__Region_5__20240628__212431", "Mild")

## ---- Merge & basic QC ----
xenium.obj <- merge(s1, y = c(s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13))
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 10)

safe_pdf("xenium_QC.pdf", 7, 5, {
  VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
})

## ---- Normalize/Reduce/Cluster ----
xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", verbose = TRUE)
xenium.obj <- RunPCA(xenium.obj, npcs = 50, features = rownames(xenium.obj))
xenium.obj <- FindNeighbors(xenium.obj, reduction = "pca", dims = 1:50)
xenium.obj <- RunUMAP(xenium.obj, dims = 1:50)
xenium.obj <- FindClusters(xenium.obj, resolution = c(0.1,0.2,0.3,0.4))

safe_pdf("SCT_snn_res.0.1_umap.pdf", 7, 5, {
  cols <- pal_disc(xenium.obj$SCT_snn_res.0.1)
  DimPlot(xenium.obj, group.by = "SCT_snn_res.0.1", cols = cols)
})

# If the fov name exists, make an image plot (example fov.13)
if ("fov.13" %in% names(xenium.obj@images)) {
  safe_pdf("fov.13_SCT_snn_res.0.1.pdf", 7, 5, {
    cols <- pal_disc(xenium.obj$SCT_snn_res.0.1)
    ImageDimPlot(xenium.obj, cols = cols, size = 0.4, fov = "fov.13", group.by = "SCT_snn_res.0.1")
  })
}

## ---- Markers ----
DefaultAssay(xenium.obj) <- "SCT"
Idents(xenium.obj) <- "SCT_snn_res.0.1"
rna.rnamarkers <- FindAllMarkers(xenium.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, "./DE_SCT_snn_res.0.1.csv", quote = FALSE, row.names = FALSE)

# Example: highlight cluster 2 on the image (guard if exists)
if (2 %in% levels(xenium.obj$SCT_snn_res.0.1)) {
  ImageDimPlot(xenium.obj, cols = "red", cells = WhichCells(xenium.obj, idents = 2))
}

## ---- Save integrated object ----
saveRDS(xenium.obj, "./xenium.obj.integrated.rds")

## ---- Example feature visuals ----
if ("fov.9" %in% names(xenium.obj@images)) {
  ImageFeaturePlot(xenium.obj, features = c("MGP"), size = 0.75, cols = c("white", "red"),
                   fov = "fov.9", max.cutoff = 3)
}

safe_pdf("dotplot_mgp_vcan_ltbp2.pdf", 4, 3.5, {
  DotPlot(xenium.obj, features = c("MGP","VCAN","LTBP2"), group.by = "disease") + RotatedAxis()
})

############################################################
# Reference mapping to a global human reference (optional)
# (Set your reference RDS path; mapping will be skipped if not found)
############################################################
ref_path <- "PATH/TO/global_reference.rds"
if (file.exists(ref_path)) {
  ref <- readRDS(ref_path)
  DefaultAssay(ref) <- "SCT"

  anchors <- FindTransferAnchors(
    reference = ref,
    query = xenium.obj,
    normalization.method = "SCT",
    reference.reduction = if ("harmony" %in% names(ref@reductions)) "harmony" else "pca",
    dims = 1:50,
    recompute.residuals = FALSE
  )

  preds_celltype <- TransferData(anchorset = anchors, refdata = ref$cell.type, dims = 1:50)
  xenium.obj <- AddMetaData(xenium.obj, metadata = preds_celltype)

  # Transfer ADT if present in ref
  if ("ADT" %in% Assays(ref)) {
    preds_adt <- TransferData(anchorset = anchors, refdata = GetAssayData(ref, assay = "ADT"), dims = 1:50)
    xenium.obj@assays[["CITE"]] <- preds_adt
  }

  safe_pdf("mapped_cell.type.disease_umap.pdf", 10, 5, {
    cols <- pal_disc(ref$cell.type)
    DimPlot(xenium.obj, group.by = "predicted.id", cols = cols, split.by = "disease")
  })

  DimPlot(xenium.obj, group.by = "predicted.id", cols = pal_disc(ref$cell.type))

  safe_pdf("fov_predicted.id.pdf", 7, 5, {
    ImageDimPlot(xenium.obj, group.by = "predicted.id",
                 cols = pal_disc(ref$cell.type), fov = "fov")
  })

  ImageFeaturePlot(xenium.obj, features = c("MFGE8"), size = 0.75, cols = c("white","red"), fov = "fov.10")

  Idents(xenium.obj) <- "predicted.id"
  if ("ModSMC" %in% levels(xenium.obj$predicted.id) && "fov.6" %in% names(xenium.obj@images)) {
    ImageDimPlot(xenium.obj, cols = "red", cells = WhichCells(xenium.obj, idents = "ModSMC"), fov = "fov.6")
  }

  write.csv(xenium.obj@meta.data, "./xenium.obj.integrated.mapped.meta.csv", quote = TRUE)
}

## ---- Save (and reload if needed) ----
# xenium.obj <- readRDS("xenium.obj.integrated.rds")
# FeaturePlot(xenium.obj, features = "OLR1")

############################################################
# Subset major compartments and write meta (if previously saved)
############################################################
# (If you already produced these sub-objects elsewhere)
if (file.exists("xenium.Endothelium.rds")) {
  xenium.Endothelium <- readRDS("xenium.Endothelium.rds")
  write.csv(xenium.Endothelium@meta.data, "./xenium.Endothelium.mapped.meta.csv", quote = TRUE)
}
if (file.exists("xenium.Myeloid.rds")) {
  xenium.Myeloid <- readRDS("xenium.Myeloid.rds")
  write.csv(xenium.Myeloid@meta.data, "./xenium.Myeloid.mapped.meta.csv", quote = TRUE)
}
if (file.exists("xenium.stroma.rds")) {
  xenium.stroma <- readRDS("xenium.stroma.rds")
  write.csv(xenium.stroma@meta.data, "./xenium.stroma.mapped.meta.csv", quote = TRUE)
}
if (file.exists("xenium.TCells.rds")) {
  xenium.TCells <- readRDS("xenium.TCells.rds")
  write.csv(xenium.TCells@meta.data, "./xenium.TCells.mapped.meta.csv", quote = TRUE)
}

############################################################
# STROMA: subset & map to SMC/Fib reference
############################################################
if ("predicted.id" %in% colnames(xenium.obj@meta.data)) {
  Idents(xenium.obj) <- "predicted.id"
}
stroma <- subset(xenium.obj, idents = c("Fibroblast2","Fibroblast1","SMCPericyte","ModSMC"))

ref_smc_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/cell_types/smc_fib/smc_fib_annotated.rds"
if (file.exists(ref_smc_path)) {
  reference <- readRDS(ref_smc_path)
  DefaultAssay(reference) <- "SCT"; DefaultAssay(stroma) <- "SCT"

  anchors <- FindTransferAnchors(
    reference = reference, query = stroma,
    normalization.method = "SCT", dims = 1:50,
    reference.assay = "SCT", reference.reduction = "pca"
  )

  stroma <- MapQuery(
    anchorset = anchors, query = stroma, reference = reference,
    refdata = list(celltype = "cell.state", predicted_ADT = "ADT"),
    reference.reduction = "pca", reduction.model = "rna.umap"
  )

  DimPlot(stroma, reduction = "ref.umap", group.by = "predicted.celltype",
          label.size = 4, cols = pal_disc(reference$cell.state), label = FALSE)

  stroma <- RunUMAP(stroma, dims = 1:50)
  RidgePlot(stroma, features = "prediction.score.max", group.by = "orig.ident")

  safe_pdf("mapped_stroma_cell.state_umap.pdf", 7, 5, {
    DimPlot(stroma, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state))
  })

  safe_pdf("fov_stroma_cell.state.pdf", 7, 5, {
    ImageDimPlot(stroma, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state), fov = "fov")
  })
}

saveRDS(stroma, "xenium.stroma.rds")

## ---- Stroma recluster & markers ----
stroma <- RunUMAP(stroma, dims = 1:50)
stroma <- FindNeighbors(stroma, reduction = "pca", dims = 1:50)
stroma <- FindClusters(stroma, resolution = c(0.1,0.2,0.3,0.4,0.5))

DimPlot(stroma, group.by = "SCT_snn_res.0.2", cols = pal_disc(stroma$SCT_snn_res.0.2))

if ("fov.10" %in% names(stroma@images)) {
  ImageDimPlot(stroma, group.by = "SCT_snn_res.0.2", cols = pal_disc(stroma$SCT_snn_res.0.2), fov = "fov.10")
}

DefaultAssay(stroma) <- "SCT"; Idents(stroma) <- "SCT_snn_res.0.2"
stroma_markers <- FindAllMarkers(stroma, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
dir.create("./stroma.cell.state", showWarnings = FALSE)
write.csv(stroma_markers, "./stroma.cell.state/DE_stroma_SCT_snn_res.0.2.csv", quote = FALSE, row.names = FALSE)

############################################################
# MYELOID: subset & map
############################################################
Myeloid <- subset(xenium.obj, idents = "Myeloid")
ref_my_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/cell_types/myeloid/final/myeloid_annotated.rds"
if (file.exists(ref_my_path)) {
  reference <- readRDS(ref_my_path)
  DefaultAssay(reference) <- "SCT"; DefaultAssay(Myeloid) <- "SCT"

  anchors <- FindTransferAnchors(reference = reference, query = Myeloid,
                                 normalization.method = "SCT", dims = 1:50,
                                 reference.assay = "SCT", reference.reduction = "pca")
  reference <- RunUMAP(reference, reduction = "harmony_rna", dims = 1:50,
                       reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
                       assay = "RNA", return.model = TRUE)
  Myeloid <- MapQuery(anchorset = anchors, query = Myeloid, reference = reference,
                      refdata = list(celltype = "cell.state", predicted_ADT = "ADT"),
                      reference.reduction = "pca", reduction.model = "rna.umap")

  DimPlot(Myeloid, reduction = "ref.umap", group.by = "predicted.celltype",
          label.size = 4, cols = pal_disc(reference$cell.state), label = FALSE)

  Myeloid <- RunUMAP(Myeloid, dims = 1:50)
  RidgePlot(Myeloid, features = "prediction.score.max", group.by = "orig.ident")

  safe_pdf("mapped_stroma_cell.state_umap.pdf", 7, 5, {
    DimPlot(Myeloid, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state))
  })

  if ("fov.13" %in% names(Myeloid@images)) {
    safe_pdf("fov.13_Myeloid_cell.state.pdf", 7, 5, {
      ImageDimPlot(Myeloid, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state), fov = "fov.13")
    })
  }

  saveRDS(Myeloid, "xenium.Myeloid.rds")
}

############################################################
# ENDOTHELIUM: subset & map
############################################################
Endothelium <- subset(xenium.obj, idents = "Endothelium")
ref_endo_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/cell_types/endothelium/endothelium.rds"
if (file.exists(ref_endo_path)) {
  reference <- readRDS(ref_endo_path)
  DefaultAssay(reference) <- "SCT"; DefaultAssay(Endothelium) <- "SCT"

  anchors <- FindTransferAnchors(reference = reference, query = Endothelium,
                                 normalization.method = "SCT", dims = 1:50,
                                 reference.assay = "SCT", reference.reduction = "pca")
  reference <- RunUMAP(reference, reduction = "harmony_rna", dims = 1:50,
                       reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
                       assay = "RNA", return.model = TRUE)
  Endothelium <- MapQuery(anchorset = anchors, query = Endothelium, reference = reference,
                          refdata = list(celltype = "cell.state", predicted_ADT = "ADT"),
                          reference.reduction = "pca", reduction.model = "rna.umap")

  DimPlot(Endothelium, reduction = "ref.umap", group.by = "predicted.celltype",
          label.size = 4, cols = pal_disc(reference$cell.state), label = FALSE)

  Endothelium <- RunUMAP(Endothelium, dims = 1:50)
  RidgePlot(Endothelium, features = "prediction.score.max", group.by = "orig.ident")

  safe_pdf("mapped_Endothelium_cell.state_umap.pdf", 7, 5, {
    DimPlot(Endothelium, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state))
  })

  if ("fov.13" %in% names(Endothelium@images)) {
    safe_pdf("fov.13_Endothelium_cell.state.pdf", 7, 5, {
      ImageDimPlot(Endothelium, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state), fov = "fov.13")
    })
  }
  saveRDS(Endothelium, "xenium.Endothelium.rds")
}

############################################################
# T/NK cells: subset & map
############################################################
TCells <- subset(xenium.obj, idents = "TCells")
ref_tnk_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/cell_types/t_nk_cells/tnkcells.rds"
if (file.exists(ref_tnk_path)) {
  reference <- readRDS(ref_tnk_path)
  DefaultAssay(reference) <- "SCT"; DefaultAssay(TCells) <- "SCT"

  anchors <- FindTransferAnchors(reference = reference, query = TCells,
                                 normalization.method = "SCT", dims = 1:50,
                                 reference.assay = "SCT", reference.reduction = "pca")
  reference <- RunUMAP(reference, reduction = "harmony_rna", dims = 1:50,
                       reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
                       assay = "RNA", return.model = TRUE)
  TCells <- MapQuery(anchorset = anchors, query = TCells, reference = reference,
                     refdata = list(celltype = "cell.state.l2", predicted_ADT = "ADT"),
                     reference.reduction = "pca", reduction.model = "rna.umap")

  DimPlot(TCells, reduction = "ref.umap", group.by = "predicted.celltype",
          label.size = 4, cols = pal_disc(reference$cell.state.l2), label = FALSE)

  TCells <- RunUMAP(TCells, dims = 1:50)
  RidgePlot(TCells, features = "prediction.score.max", group.by = "orig.ident")

  safe_pdf("mapped_TCells_cell.state_umap.pdf", 7, 5, {
    DimPlot(TCells, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state.l2))
  })

  if ("fov.13" %in% names(TCells@images)) {
    safe_pdf("fov.13_TCells_cell.state.pdf", 7, 5, {
      ImageDimPlot(TCells, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state.l2), fov = "fov.13")
    })
  }

  saveRDS(TCells, "xenium.TCells.rds")
}

############################################################
# Stroma inspection & ridge summary
############################################################
stroma <- readRDS("xenium.stroma.rds")
if (exists("reference") && "cell.state" %in% colnames(reference@meta.data)) {
  DimPlot(stroma, group.by = "predicted.celltype", cols = pal_disc(reference$cell.state))
  safe_pdf("mapping_ridge.pdf", 6, 4, {
    RidgePlot(stroma, features = "predicted.celltype.score", group.by = "predicted.celltype") +
      scale_fill_manual(values = pal_disc(reference$cell.state))
  })
}

############################################################
# Key spatial features across fovs
############################################################
genesA <- c("MGP","VCAN","LTBP2","TNFRSF11B")
for (g in genesA) {
  for (f in names(xenium.obj@images)) {
    try({
      ImageFeaturePlot(xenium.obj, features = g, size = 0.75, cols = c("white","red"),
                       fov = f, max.cutoff = ifelse(g %in% c("LTBP2","TNFRSF11B"), 2, 3))
    }, silent = TRUE)
  }
}

# Optional: FAP-like module intersection (against current object genes)
FAP_module <- c("MGP","LUM","VCAN","F2R","OMD","LTBP2","FAP","FN1","COL1A2","CFH","PLXDC2",
                "THBS2","ITGBL1","PRELP","FBLN2","PDFRA","DCN","MMP2","PDGFD,","FXYD5",
                "POSTN","LAMA2","MEG3","COL3A1","TMSB10","CRTAC1","PCOLCE2","COL8A1")
x <- intersect(FAP_module, rownames(xenium.obj))
# (You can AddModuleScore on SCT if desired)
# DefaultAssay(xenium.obj) <- "SCT"
# xenium.obj <- AddModuleScore(xenium.obj, features = list(x), name = "FAPModule")
