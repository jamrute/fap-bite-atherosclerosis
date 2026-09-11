################################################################################
# Human CITE-seq: WNN (RNA + protein) integration and global cell-type annotation
#
# Paper : Amrute et al., Science (2026) | doi:10.1126/science.adx1736
# Part  : Human coronary artery CITE-seq atlas
#
# Purpose
#   Builds a weighted-nearest-neighbour graph from RNA and ADT embeddings,
#   annotates clusters (SCT_snn_res.0.2) into major cell types, finds RNA and
#   protein markers per cell type, plots marker heatmaps and cell-type
#   composition by sample/clinical metadata, and exports objects for Scanpy
#   and for Tangram deconvolution of the Visium data.
#
# Inputs
#   RNA_clustered_postContamination.rds  (normalized & clustered object with
#     SCT and ADT assays, pca/apca, rna.umap/adt.umap; see Notes)
#
# Outputs
#   DE_SCT_cell.type.csv, DE_ADT_cell.type.csv  (cell-type markers)
#   heatmap_GEX.png, heatmap_ADT.png
#   annotated.h5Seurat / annotated.h5ad
#   annotated_raw.h5Seurat / annotated_raw.h5ad + meta_raw.csv  (Tangram reference)
#
# Run order
#   Upstream  : 01_merge_qc_doublet_filter.R (+ normalization step, see Notes)
#   Downstream: 04, 05, 08, spatial/visium/01_tangram_deconvolution.ipynb
#
# Notes
#   Interactive analysis script: run section by section (e.g. in RStudio).
#   The step between 01 and this script (SCTransform, ADT normalization,
#   PCA/APCA, contamination clean-up, initial clustering) is not included.
#   Downstream scripts read the annotated object as `integrated_annotated.rds`.
################################################################################

## ---- Shared helpers (R/utils.R) ----
repo_dir <- "."  # EDIT: path to the root of this repository
source(file.path(repo_dir, "R", "utils.R"))

## ---- Libraries (unique + used here) ----
library(Seurat)       # core single-cell workflow (WNN, UMAP, clustering, markers)
library(SeuratDisk)   # Save/Convert Seurat <-> h5ad
library(dplyr)        # data wrangling (joins for metadata)
library(ggplot2)      # plots
library(patchwork)    # plot composition
library(harmony)      # not used directly in this script
library(ArchR)        # for paletteContinuous() used in heatmaps/FeaturePlot
library(ggsci)        # additional colour palettes
library(scProportionTest)  # proportion tests
library(cowplot)      # plotting helpers
library(Nebulosa)     # plot_density()

## ============================================================
## Load clustered object
## ============================================================
cor <- readRDS("./RNA_clustered_postContamination.rds")

## Quick cluster UMAP
DimPlot(cor, reduction = 'rna.umap', group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)

## ============================================================
## WNN graph: RNA (PCA) + ADT (APCA)
## ============================================================
cor <- FindMultiModalNeighbors(
  object = cor,
  reduction.list = list("pca", "apca"),
  weighted.nn.name = "sct.dsb_wnn",
  knn.graph.name = "sct.dsb_knn",
  modality.weight.name = "sct.dsb_weight",
  snn.graph.name = "sct.dsb_snn",
  dims.list = list(1:80, 1:30),
  verbose = TRUE
)

## UMAP on WNN graph and clustering on WNN SNN
cor <- RunUMAP(
  cor, nn.name = "sct.dsb_wnn",
  reduction.name = "sct.dsb_wnn_umap",
  reduction.key = "sct.dsb_wnnUMAP_", verbose = TRUE,
  return.model = TRUE
)

# Optional: cluster directly on the WNN graph. The annotation below uses the
# RNA-based clustering (SCT_snn_res.0.2) already stored in the input object.
# cor <- FindClusters(cor, graph.name = "sct.dsb_snn", algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = TRUE)

## Compare embeddings colored by the same clustering
DimPlot(cor, reduction = 'rna.umap',          group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)
DimPlot(cor, reduction = 'adt.umap',          group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)
DimPlot(cor, reduction = 'sct.dsb_wnn_umap',  group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)

## ============================================================
## Annotate global clusters (SCT_snn_res.0.2) -> cell.type
## ============================================================
cell_type_labels <- c(
  "0" = "TCells", "1" = "Myeloid", "2" = "SMCPericyte", "3" = "Endothelium",
  "4" = "TCells", "5" = "Fibroblast1", "6" = "BCells", "7" = "Fibroblast2",
  "8" = "TCells", "9" = "ModSMC", "10" = "Mast", "11" = "Glia",
  "12" = "Lymphatic", "13" = "PlasmaCells", "14" = "Proliferating", "15" = "pDC",
  "16" = "Myeloid"
)
cor$cell.type <- annotate_clusters(cor$SCT_snn_res.0.2, cell_type_labels)

cor$cell.type <- factor(
  cor$cell.type,
  levels = c("SMCPericyte","ModSMC","Fibroblast1","Fibroblast2","Endothelium",
             "Myeloid","pDC","Mast","BCells","PlasmaCells","TCells",
             "Proliferating","Lymphatic","Glia")
)

## Cell type UMAPs across modalities
DimPlot(cor, reduction = 'rna.umap', group.by = 'cell.type', label.size = 4,
        cols = paletteDiscrete(unique(cor$cell.type), set = "stallion"), label = FALSE)
DimPlot(cor, reduction = 'adt.umap', group.by = 'cell.type', label.size = 4,
        cols = paletteDiscrete(unique(cor$cell.type), set = "stallion"), label = FALSE)
DimPlot(cor, reduction = 'sct.dsb_wnn_umap', group.by = 'cell.type', label.size = 4,
        cols = paletteDiscrete(unique(cor$cell.type), set = "stallion"), label = FALSE)

## ============================================================
## DGE per cell.type (SCT then ADT)
## ============================================================
DefaultAssay(cor) <- 'SCT'
Idents(cor) <- "cell.type"
rna.rnamarkers <- FindAllMarkers(cor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./DE_SCT_cell.type.csv", quote = FALSE)

DefaultAssay(cor) <- 'ADT'
Idents(cor) <- "cell.type"
rna.rnamarkers <- FindAllMarkers(cor, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./DE_ADT_cell.type.csv", quote = FALSE)

## ============================================================
## Marker heatmaps (GEX and ADT)
## ============================================================
# GEX markers
genes <- c("MYH11","CCL19","FBLN1","PTGDS","VWF","C1QA","GZMB","KIT","CD79A","JCHAIN","CD3E","MKI67","CCL21","NRXN1")
DoHeatmap(cor, features = genes, assay = "SCT",
          group.colors = paletteDiscrete(unique(cor$cell.type), set = "stallion"),
          size = 2, angle = 90) + NoLegend() +
  scale_fill_gradientn(colours = paletteContinuous("solarExtra"))
ggsave(filename = "heatmap_GEX.png")

# ADT markers
adts <- c("ITGB1.1","FAP.1","THY1.1","MME.1","PECAM1.1","CSF2RA.1","IL3RA.1",
          "KIT.1","CD19.1","CD27.1","CD3D.1","CD2.1","MRC1.1","NCAM1.1")
DoHeatmap(cor, features = adts, assay = "ADT",
          group.colors = paletteDiscrete(unique(cor$cell.type), set = "stallion"),
          size = 2, angle = 90) + NoLegend() +
  scale_fill_gradientn(colours = paletteContinuous("blueYellow"))
ggsave(filename = "heatmap_ADT.png")

## Example density plot
DefaultAssay(cor) <- "SCT"
plot_density(cor, features = "ITGA9", reduction = 'rna.umap')

## ============================================================
## Cell-type composition by sample and clinical metadata
## ============================================================
plot_composition(cor@meta.data, x = "sampleID", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

plot_composition(cor@meta.data, x = "AgeRange", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

plot_composition(cor@meta.data, x = "Disease", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

plot_composition(cor@meta.data, x = "Stent", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

plot_composition(cor@meta.data, x = "HF", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

plot_composition(cor@meta.data, x = "Sex", fill = "cell.type",
                 colors = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion")))

## ============================================================
## Export h5Seurat/h5ad for Scanpy; keep SCT model median_umi
## ============================================================
cor$cell.type <- as.character(cor$cell.type) 
slot(cor$SCT@SCTModel.list[[1]], 'median_umi') <- median(cor$SCT@SCTModel.list[[1]]@cell.attributes$umi)
SaveH5Seurat(cor, filename = "./annotated.h5Seurat")
Convert("./annotated.h5Seurat", dest = "h5ad")

## ============================================================
## Export raw RNA counts only (Tangram reference for spatial/visium/01)
## ============================================================
DefaultAssay(cor) <- "RNA"
cor@assays[["SCT"]] <- NULL
cor@assays[["ADT"]] <- NULL
SaveH5Seurat(cor, filename = "./annotated_raw.h5Seurat")
Convert("./annotated_raw.h5Seurat", dest = "h5ad")
write.csv(cor@meta.data, file = "./meta_raw.csv", quote = TRUE)

