###############################################################
# CITE-seq: WNN integration, annotation, markers, got global
# CITE-seq object
###############################################################

## ---- Libraries (unique + used here) ----
library(Seurat)       # core single-cell workflow (WNN, UMAP, clustering, markers)
library(SeuratDisk)   # Save/Convert Seurat <-> h5ad
library(dplyr)        # data wrangling (joins for metadata)
library(ggplot2)      # plots
library(patchwork)    # plot composition
library(harmony)      # (kept from original list; not used below but safe to keep)
library(ArchR)        # for paletteContinuous() used in heatmaps/FeaturePlot
library(ggsci)        # palettes if you switch from paletteDiscrete()
library(scProportionTest)  # proportion tests
library(cowplot)      # plotting helpers (used near the end)

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

# NOTE: The original code used `sample <- FindClusters(sample, ...)`.
# If your object is `cor`, run the following line instead (left as a comment to avoid changing code):
# cor <- FindClusters(cor, graph.name = "sct.dsb_snn", algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = TRUE)

## Compare embeddings colored by the same clustering
DimPlot(cor, reduction = 'rna.umap',          group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)
DimPlot(cor, reduction = 'adt.umap',          group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)
DimPlot(cor, reduction = 'sct.dsb_wnn_umap',  group.by = 'SCT_snn_res.0.2', label.size = 4, label = TRUE)

## ============================================================
## Annotate global clusters -> cell.type (kept as in your code)
## ============================================================
fun <- function(x) {
  if (x == "0") {"TCells"} 
  else if (x == "1") {"Myeloid"}
  else if (x == "2") {"SMCPericyte"}
  else if (x == "3") {"Endothelium"}
  else if (x == "4") {"TCells"}
  else if (x == "5") {"Fibroblast1"}
  else if (x == "6") {"BCells"}
  else if (x == "7") {"Fibroblast2"}
  else if (x == "8") {"TCells"}
  else if (x == "9") {"ModSMC"}
  else if (x == "10") {"Mast"}
  else if (x == "11") {"Glia"}
  else if (x == "12") {"Lymphatic"}
  else if (x == "13") {"PlasmaCells"}
  else if (x == "14") {"Proliferating"}
  else if (x == "15") {"pDC"}
  else if (x == "16") {"Myeloid"}
}
cor$cell.type <- mapply(fun, cor$SCT_snn_res.0.2)

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
## Composition plots by metadata (unchanged from your code)
## ============================================================
ggplot(cor@meta.data, aes(x = sampleID, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(cor@meta.data, aes(x = AgeRange, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(cor@meta.data, aes(x = Disease, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(cor@meta.data, aes(x = Stent, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(cor@meta.data, aes(x = HF, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(cor@meta.data, aes(x = Sex, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(cor$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

## ============================================================
## Export h5Seurat/h5ad for Scanpy; keep SCT model median_umi
## ============================================================
cor$cell.type <- as.character(cor$cell.type) 
slot(cor$SCT@SCTModel.list[[1]], 'median_umi') <- median(cor$SCT@SCTModel.list[[1]]@cell.attributes$umi)
SaveH5Seurat(cor, filename = "./annotated.h5Seurat")
Convert("./annotated.h5Seurat", dest = "h5ad")

## ============================================================
## Export "raw-ish" RNA only for Tangram
## ============================================================
DefaultAssay(cor) <- "RNA"
cor@assays[["SCT"]] <- NULL
cor@assays[["ADT"]] <- NULL
SaveH5Seurat(cor, filename = "./annotated_raw.h5Seurat")
Convert("./annotated_raw.h5Seurat", dest = "h5ad")
write.csv(cor@meta.data, file = "./meta_raw.csv", quote = TRUE)

