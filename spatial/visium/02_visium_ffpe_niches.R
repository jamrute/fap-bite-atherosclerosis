############################################################
# Visium FFPE (multi-sample) spatial analysis + Tangram composition
# - Libraries deduplicated & ordered
# - Adds missing helpers (e.g., cor.mtest)
# - Preserves your exact workflow, with small safety tweaks
############################################################

## ---- Libraries (unique) ----
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ArchR)          # paletteDiscrete / paletteContinuous
library(Matrix)
library(data.table)
library(gt)
library(igraph)
library(RColorBrewer)
library(compositions)
library(harmony)
library(tidyverse)
library(clustree)
library(uwot)
library(scran)
library(cluster)
library(pheatmap)
library(corrplot)
library(remotes)
# remotes::install_github("carmonalab/GeneNMF")
library(GeneNMF)
library(UCell)
library(RcppML)
library(viridis)
library(scCustomize)
library(msigdbr)
library(fgsea)
library(UpSetR)

## ---- Helper: corrplot p-value test ----
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

############################################################
# 1) Select cells in T1099R (manual exclusion), preprocess all slides
############################################################

# --- T1099R (interactive selection & keep) ---
s6 <- Load10X_Spatial(data.dir = "../data/Visium_T1099R", slice = "Visium_T1099R")
s6$sample <- "Visium_T1099R"; s6$category <- "Mild"

spot_composition <- read.csv("../tangram/output/Visium_T1099R/Visium_T1099R.csv")
spot_composition2 <- spot_composition[, -1]; rownames(spot_composition2) <- spot_composition[, 1]
s6$cell.type.max <- apply(spot_composition2, 1, function(x) paste0(names(spot_composition2)[x == max(x)], collapse = "_"))
spot_composition2 <- t(spot_composition2)
s6@assays[["SpotComposition"]] <- CreateAssayObject(spot_composition2)
s6@assays$SpotComposition@key <- "spotcomposition_"

s6 <- subset(s6, subset = nCount_Spatial > 0)
s6 <- SCTransform(s6, assay = "Spatial", verbose = FALSE, vars.to.regress = "nCount_Spatial")
DefaultAssay(s6) <- "SCT"
s6 <- RunPCA(s6, assay = "SCT", verbose = FALSE)
s6 <- FindNeighbors(s6, reduction = "pca", dims = 1:30)
s6 <- FindClusters(s6, verbose = FALSE, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))
s6 <- RunUMAP(s6, reduction = "pca", dims = 1:30)

DimPlot(s6, reduction = "umap", group.by = "SCT_snn_res.0.5",
        cols = paletteDiscrete(unique(s6$SCT_snn_res.0.5), set = "stallion"))

Idents(s6) <- "SCT_snn_res.0.5"
SpatialDimPlot(s6, cols = paletteDiscrete(unique(s6$SCT_snn_res.0.5), set = "stallion"))

# Manual point selection (click on the plot window created by ggplot)
df <- s6@images[["Visium_T1099R"]]@coordinates
df$cluster <- as.character(s6$SCT_snn_res.0.5)
x1 <- CellSelector(df %>% ggplot(aes(x = imagerow, y = imagecol, color = cluster)) + geom_point())

s6 <- subset(s6, cells = x1, invert = TRUE)
s6_keep <- s6

# --- Helper to load & annotate a slide quickly ---
load_visium <- function(datadir, slice, sample_name, category, tangram_csv) {
  obj <- Load10X_Spatial(data.dir = datadir, slice = slice)
  obj$sample <- sample_name
  obj$category <- category
  spot_comp <- read.csv(tangram_csv)
  sc2 <- spot_comp[, -1]; rownames(sc2) <- spot_comp[, 1]
  obj$cell.type.max <- apply(sc2, 1, function(x) paste0(names(sc2)[x == max(x)], collapse = "_"))
  sc2 <- t(sc2)
  obj@assays[["SpotComposition"]] <- CreateAssayObject(sc2)
  obj@assays$SpotComposition@key <- "spotcomposition_"
  obj <- subset(obj, subset = nCount_Spatial > 0)
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE, vars.to.regress = "nCount_Spatial")
  obj
}

# --- Load remaining slides ---
s1  <- load_visium("../data/Visium_D1118_LAD", "Visium_D1118_LAD", "Visium_D1118_LAD", "Severe",   "../tangram/output/Visium_D1118_LAD/Visium_D1118_LAD.csv")
s2  <- load_visium("../data/Visium_LCA",       "Visium_LCA",       "Visium_LCA",       "Severe",   "../tangram/output/Visium_LCA/Visium_LCA.csv")
s3  <- load_visium("../data/Visium_T1071",     "Visium_T1071",     "Visium_T1071",     "Severe",   "../tangram/output/Visium_T1071/Visium_T1071.csv")
s4  <- load_visium("../data/Visium_T1077",     "Visium_T1077",     "Visium_T1077",     "Moderate", "../tangram/output/Visium_T1077/Visium_T1077.csv")
s5  <- load_visium("../data/Visium_T1079",     "Visium_T1079",     "Visium_T1079",     "Moderate", "../tangram/output/Visium_T1079/Visium_T1079.csv")
s7  <- load_visium("../data/Visium_T1111",     "Visium_T1111",     "Visium_T1111",     "Severe",   "../tangram/output/Visium_T1111/Visium_T1111.csv")
s8  <- load_visium("../data/Visium_T1153R",    "Visium_T1153R",    "Visium_T1153R",    "Moderate", "../tangram/output/Visium_T1153R/Visium_T1153R.csv")
s9  <- load_visium("../data/Visium_T1161_RCA", "Visium_T1161_RCA", "Visium_T1161_RCA", "Moderate", "../tangram/output/Visium_T1161_RCA/Visium_T1161_RCA.csv")
s10 <- load_visium("../data/Visium_T1163",     "Visium_T1163",     "Visium_T1163",     "Moderate", "../tangram/output/Visium_T1163/Visium_T1163.csv")
s11 <- load_visium("../data/Visium_T1167_LAD", "Visium_T1167_LAD", "Visium_T1167_LAD", "Severe",   "../tangram/output/Visium_T1167_LAD/Visium_T1167_LAD.csv")
s12 <- load_visium("../data/Visium_T1170_LAD", "Visium_T1170_LAD", "Visium_T1170_LAD", "Mild",     "../tangram/output/Visium_T1170_LAD/Visium_T1170_LAD.csv")
s13 <- load_visium("../data/Visium_T1171_LM",  "Visium_T1171_LM",  "Visium_T1171_LM",  "Mild",     "../tangram/output/Visium_T1171_LM/Visium_T1171_LM.csv")
s14 <- load_visium("../data/Visium_T1183_RCA", "Visium_T1183_RCA", "Visium_T1183_RCA", "Moderate", "../tangram/output/Visium_T1183_RCA/Visium_T1183_RCA.csv")
s15 <- load_visium("../data/Visium_T1186_RCA", "Visium_T1186_RCA", "Visium_T1186_RCA", "Mild",     "../tangram/output/Visium_T1186_RCA/Visium_T1186_RCA.csv")
s16 <- load_visium("../data/Visium_T1190_RCA", "Visium_T1190_RCA", "Visium_T1190_RCA", "Moderate", "../tangram/output/Visium_T1190_RCA/Visium_T1190_RCA.csv")

# Merge all
merged <- merge(s1, y = c(s2, s3, s4, s5, s6_keep, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16))
saveRDS(merged, "merged.rds")

# Or reload integrated object
# merged <- readRDS("cad_spatial_integrated.rds")

############################################################
# 2) RNA (SCT) PCA/Neighbors/Clusters/UMAP + spatial plots
############################################################
DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- c(VariableFeatures(s1), VariableFeatures(s2), VariableFeatures(s3),
                              VariableFeatures(s4), VariableFeatures(s5), VariableFeatures(s6_keep),
                              VariableFeatures(s7), VariableFeatures(s8), VariableFeatures(s9),
                              VariableFeatures(s10), VariableFeatures(s11), VariableFeatures(s12),
                              VariableFeatures(s13), VariableFeatures(s14), VariableFeatures(s15),
                              VariableFeatures(s16))
merged <- RunPCA(merged, assay = "SCT", verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:50)
merged <- FindClusters(merged, verbose = FALSE, resolution = c(0.1,0.2,0.3,0.4,0.5))
merged <- RunUMAP(merged, reduction = "pca", dims = 1:50)

DimPlot(merged, reduction = "umap", group.by = "SCT_snn_res.0.3",
        cols = paletteDiscrete(unique(merged$SCT_snn_res.0.3), set = "stallion"), split.by = "category")

DimPlot(merged, reduction = "umap", group.by = "SCT_snn_res.0.3",
        cols = paletteDiscrete(unique(merged$SCT_snn_res.0.3), set = "stallion"), ncol = 2)

Idents(merged) <- "SCT_snn_res.0.3"

# One PNG per slide
for (im in unique(merged$sample)) {
  try({
    SpatialDimPlot(merged, images = im, cols = paletteDiscrete(unique(merged$SCT_snn_res.0.3), set = "stallion"))
    ggsave(paste0(im, "-SCT_snn_res.0.3.png"), dpi = 300, width = 6, height = 5)
  }, silent = TRUE)
}

# DE (SCT)
merged <- PrepSCTFindMarkers(merged, assay = "SCT", verbose = TRUE)
DefaultAssay(merged) <- "SCT"
Idents(merged) <- "SCT_snn_res.0.3"
rna.rnamarkers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(rna.rnamarkers, file = "./DE_SCT_snn_res.0.3.csv", quote = FALSE)

############################################################
# 3) Spatially Variable Genes (Moran's I)
############################################################
merged <- FindSpatiallyVariableFeatures(merged, assay = "SCT",
                                        features = VariableFeatures(merged)[1:3000],
                                        selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(merged, selection.method = "moransi"), 6)

for (im in unique(merged$sample)) {
  try({
    SpatialFeaturePlot(merged, features = top.features, ncol = 3, images = im)
    ggsave(paste0(im, "-Top6_SVG.png"), dpi = 300, width = 7, height = 6)
  }, silent = TRUE)
}

############################################################
# 4) PROGENy scores (NFkB / TGFb) + spatial
############################################################
Idents(merged) <- "SCT_snn_res.0.3"
CellsClusters <- data.frame(Cell = names(Idents(merged)),
                            CellType = as.character(Idents(merged)),
                            stringsAsFactors = FALSE)

merged <- progeny(merged, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
merged <- Seurat::ScaleData(merged, assay = "progeny")

DefaultAssay(merged) <- "progeny"
FeaturePlot(merged, features = "NFkB", pt.size = .5, reduction = 'umap') +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"), oob = scales::squish, limits = c(15, 50)) + NoLegend()
FeaturePlot(merged, features = "TGFb", pt.size = .5, reduction = 'umap') +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"), oob = scales::squish, limits = c(25, 50)) + NoLegend()

images <- unique(merged$sample)
for (im in images) {
  try({
    SpatialFeaturePlot(merged, features = "NFkB", images = im) +
      scale_fill_gradientn(colors = paletteContinuous("blueYellow"), oob = scales::squish, limits = c(25, 50))
    ggsave(paste0("./progeny/", im, "-NFkB-Top6_SVG.png"), dpi = 300, width = 6, height = 5)

    SpatialFeaturePlot(merged, features = "TGFb", images = im) +
      scale_fill_gradientn(colors = paletteContinuous("blueYellow"), oob = scales::squish, limits = c(25, 50))
    ggsave(paste0("./progeny/", im, "-TGFb-Top6_SVG.png"), dpi = 300, width = 6, height = 5)
  }, silent = TRUE)
}

############################################################
# 5) Whole-dataset correlation (Tangram SpotComposition)
############################################################
df <- t(as.data.frame(merged@assays$SpotComposition@counts))
df <- df[, c('Fibroblast1','Myeloid','Endothelium','SMCPericyte','ModSMC','Fibroblast2','TCells','BCells')]
corr_mat <- cor(df)
pheatmap(corr_mat, scale = "none", cluster_rows = FALSE, fontsize_row = 8, fontsize_col = 8,
         cluster_cols = FALSE, legend = TRUE, cellwidth = 16, cellheight = 16, border_color = NA)

diag(corr_mat) <- 0
testRes <- cor.mtest(df, conf.level = 0.95)
corrplot(corr_mat, title = "\n\n\n", is.corr = FALSE, order = "hclust",
         diag = FALSE, method = 'color', pch.col = 'grey20',
         addrect = 3, rect.col = 'black', rect.lwd = 3)

avg <- AverageExpression(merged, assay = "SpotComposition", group.by = "SCT_snn_res.0.3")
pheatmap(as.matrix(avg$SpotComposition), scale = "row", cluster_rows = TRUE, fontsize_row = 8, fontsize_col = 8,
         cluster_cols = TRUE, legend = TRUE, cellwidth = 16, cellheight = 16, border_color = NA,
         col = viridis(240), breaks = seq(0, 1, by = 1/240))

############################################################
# 6) NMF on SCT, meta-program discovery, and MP signatures
############################################################
# (Optionally) re-SCT for uniformity
merged <- SCTransform(merged, assay = "Spatial", verbose = FALSE, vars.to.regress = "nCount_Spatial")
seu <- merged

ndim <- 30
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- runNMF(seu, k = ndim, assay = "SCT")
seu <- RunUMAP(seu, reduction = "NMF", dims = 1:ndim, reduction.name = "NMF_UMAP", reduction.key = "nmfUMAP_")

DimPlot(seu, reduction = "umap", group.by = "SCT_snn_res.0.3",
        cols = paletteDiscrete(unique(seu$SCT_snn_res.0.3), set = "stallion"), ncol = 2)

seu <- FindNeighbors(seu, reduction = "NMF", dims = 1:30, graph.name = "NMF_NN")
seu <- FindClusters(seu, verbose = FALSE, resolution = c(0.1,0.2,0.3,0.4,0.5), graph.name = "NMF_NN")
seu <- RunUMAP(seu, reduction = "NMF", dims = 1:30)

DimPlot(seu, reduction = "NMF_UMAP", group.by = "NMF_NN_res.0.5",
        cols = paletteDiscrete(unique(seu$NMF_NN_res.0.5), set = "stallion"), ncol = 2)

# Marker discovery
DefaultAssay(seu) <- "SCT"
Idents(seu) <- "SCT_snn_res.0.3"
rna.rnamarkers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(rna.rnamarkers, file = "./DE_NMF_NN_res.0.5.csv", quote = FALSE)

# Multi-NMF meta-programs
seu.list <- SplitObject(seu, split.by = "sample")
geneNMF.programs <- multiNMF(seu.list, assay = "SCT", slot = "data", k = 4:9, nfeatures = 3000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs, nMP = 10, weight.explained = 0.7, max.genes = 100)

plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)

# Collate MP genes into a single data.frame
bindMP <- function(mp_list, tag) {
  df <- data.frame(mp_list, row.names = NULL)
  colnames(df) <- "gene"; df$MP <- tag; df
}
MP_genes <- bind_rows(
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP1 , "MP1"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP2 , "MP2"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP3 , "MP3"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP4 , "MP4"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP5 , "MP5"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP6 , "MP6"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP7 , "MP7"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP8 , "MP8"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP9 , "MP9"),
  bindMP(geneNMF.metaprograms$metaprograms.genes$MP10, "MP10")
)

# GWAS overlap
GWAS <- read_delim("./CAD_GWAS_genes.csv", ",", escape_double = FALSE, trim_ws = TRUE)
gwas_genes <- GWAS$gene
intersecting_genes <- intersect(MP_genes$gene, gwas_genes)
DGE_subset <- filter(MP_genes, gene %in% intersecting_genes)
cad_gwas_overlap <- data.frame(sort(table(DGE_subset$MP), decreasing = TRUE))

pdf("./cad_gwas_overlap_MP_VisiumFFPE.pdf", useDingbats = FALSE, width = 3.5, height = 2.5)
cad_gwas_overlap$Var1 <- factor(cad_gwas_overlap$Var1,
                                levels = c("MP10","MP8","MP7","MP6","MP3","MP1","MP2","MP5","MP4"))
ggplot(cad_gwas_overlap, aes(x = Var1, y = Freq)) +
  geom_segment(aes(xend = Var1, yend = 0.5)) +
  geom_point(size = 2, color = "black") +
  coord_flip() + theme_bw() + xlab("")
dev.off()

# MP UCell scores and custom embedding
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(seu, features = mp.genes, assay = "SCT", name = "")

matrix <- seu@meta.data[, names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_", seq_len(ncol(dimred)))

seu@reductions[["MPsignatures"]] <- new("DimReduc",
  cell.embeddings = dimred, assay.used = "Spatial", key = "MP_", global = FALSE)

set.seed(123)
seu <- RunUMAP(seu, reduction = "MPsignatures",
               dims = 1:length(seu@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")

seu <- FindNeighbors(seu, reduction = "MPsignatures",
                     dims = 1:length(seu@reductions[["MPsignatures"]]),
                     graph.name = "MPsignatures_NN")
seu <- FindClusters(seu, verbose = FALSE, resolution = c(0.1,0.2,0.3,0.4,0.5),
                    graph.name = "MPsignatures_NN")

DimPlot(seu, reduction = "umap_MP", group.by = "MPsignatures_NN_res.0.5",
        cols = paletteDiscrete(unique(seu$MPsignatures_NN_res.0.5), set = "stallion"), ncol = 2)

# MP features on embeddings and spatial
FeaturePlot(seu, features = "MP4", reduction = "umap_MP") &
  scale_color_viridis(option = "B") &
  theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

FeaturePlot(seu, features = names(mp.genes), reduction = "umap_MP", ncol = 5) &
  scale_color_viridis(option = "B") &
  theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

FeaturePlot(seu, features = names(mp.genes), reduction = "NMF_UMAP", ncol = 5) &
  scale_color_viridis(option = "B") &
  theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

for (mp in paste0("MP", 1:10)) {
  try({
    SpatialFeaturePlot(seu, features = mp, images = "Visium_T1071") + scale_fill_gradientn(colors = viridis::inferno(256))
  }, silent = TRUE)
}

# MP by cell type / condition / cluster
MP_celltype <- seu@meta.data[, c("cell.type.max", paste0("MP", 1:10))]
MP_celltype_avg <- MP_celltype %>% group_by(cell.type.max) %>% summarise(across(everything(), mean))
MP_celltype_avg2 <- as.data.frame(MP_celltype_avg[, -1]); rownames(MP_celltype_avg2) <- MP_celltype_avg$cell.type.max
MP_celltype_avg2 <- MP_celltype_avg2[c("BCells","Endothelium","Fibroblast1","Fibroblast2","ModSMC","Myeloid","SMCPericyte","TCells"),]
pdf("./MP_celltype.pdf", useDingbats = FALSE, width = 3.5, height = 2.5)
pheatmap(MP_celltype_avg2, scale = "column", col = paletteContinuous("solarExtra"),
         cluster_rows = TRUE, fontsize_row = 10, fontsize_col = 10, cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

MP_cond <- seu@meta.data[, c("category", paste0("MP", 1:10))]
MP_cond_avg <- MP_cond %>% group_by(category) %>% summarise(across(everything(), mean))
MP_cond_avg2 <- as.data.frame(MP_cond_avg[, -1]); rownames(MP_cond_avg2) <- MP_cond_avg$category
pdf("./MP_severity.pdf", useDingbats = FALSE, width = 2.5, height = 1.3)
pheatmap(MP_cond_avg2, scale = "column", col = paletteContinuous("solarExtra"),
         cluster_rows = TRUE, fontsize_row = 5, fontsize_col = 5, cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

MP_cl <- seu@meta.data[, c("MPsignatures_NN_res.0.5", paste0("MP", 1:10))]
MP_cl_avg <- MP_cl %>% group_by(MPsignatures_NN_res.0.5) %>% summarise(across(everything(), mean))
MP_cl_avg2 <- as.data.frame(MP_cl_avg[, -1]); rownames(MP_cl_avg2) <- MP_cl_avg$MPsignatures_NN_res.0.5
pdf("./MP_NMF_NN_res.0.5_2.pdf", useDingbats = FALSE, width = 2.5, height = 2.5)
pheatmap(MP_cl_avg2, scale = "column", col = viridis(10), cluster_rows = TRUE, fontsize_row = 5,
         fontsize_col = 5, cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

saveRDS(seu, "integrated.rds")
saveRDS(geneNMF.metaprograms, "geneNMF.metaprograms.rds")

############################################################
# 7) FMC/CMC module (FAP-like) + MP overlays
############################################################
rna.rnamarkers <- read.csv("/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/cell_types/smc_fib/DE_RNA_cell.state.csv")
rna.rnamarkers <- rna.rnamarkers %>% filter(cluster == "CMC", avg_log2FC > 0.58)

# (Use seu gene names for intersection; user 'xenium.obj' replaced by 'seu')
x <- intersect(c("MGP","LUM","VCAN","F2R","OMD","LTBP2","FAP","FN1","COL1A2","CFH","PLXDC2",
                 "THBS2","ITGBL1","PRELP","FBLN2","PDFRA","DCN","MMP2","PDGFD,","FXYD5",
                 "POSTN","LAMA2","MEG3","COL3A1","TMSB10","CRTAC1","PCOLCE2","COL8A1"),
               rownames(seu))

DefaultAssay(seu) <- "SCT"
seu <- AddModuleScore(seu, features = list(x), assay = "SCT", name = "CMCz")

SpatialFeaturePlot(seu, features = "MP7", images = "Visium_T1071") +
  scale_fill_gradientn(colors = viridis::inferno(256))

saveRDS(seu, "integrated.rds")

Idents(seu) <- "MPsignatures_NN_res.0.5"
for (im in unique(seu$sample)) {
  try({
    SpatialDimPlot(seu, images = im, cols = paletteDiscrete(unique(seu$MPsignatures_NN_res.0.5), set = "stallion"))
    ggsave(paste0(im, "-MPsignatures_NN_res.0.5.png"), dpi = 300, width = 6, height = 5)
  }, silent = TRUE)
}

# Rebuild MP_genes later if needed
# geneNMF.metaprograms <- readRDS("geneNMF.metaprograms.rds")
write.csv(MP_genes, file = "./MP_genes.csv", quote = FALSE)

DefaultAssay(seu) <- "SCT"
pdf("./MP_dotplot.pdf", useDingbats = FALSE, width = 16, height = 4)
DotPlot(seu, features = c("MP1","MP2","MP3","MP4"), group.by = "MPsignatures_NN_res.0.5",
        col.min = 0, cols = c("lightgrey","red")) + RotatedAxis()
dev.off()

DefaultAssay(seu) <- "SCT"
Idents(seu) <- "MPsignatures_NN_res.0.5"
rna.rnamarkers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(rna.rnamarkers, file = "./DE_MPsignatures_NN_res.0.5.csv", quote = FALSE)

rna.rnamarkers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10
pdf("./MP_dotplot.pdf", useDingbats = FALSE, width = 16, height = 4)
DotPlot(seu, features = unique(top10$gene), group.by = "MPsignatures_NN_res.0.5",
        col.min = 0, cols = c("lightgrey","red")) + RotatedAxis()
dev.off()

# FAP module across MP signatures
x <- intersect(c("MGP","LUM","VCAN","F2R","OMD","LTBP2","FAP","FN1","COL1A2","CFH","PLXDC2",
                 "THBS2","ITGBL1","PRELP","FBLN2","PDFRA","DCN","MMP2","PDGFD,","FXYD5",
                 "POSTN","LAMA2","MEG3","COL3A1","TMSB10","CRTAC1","PCOLCE2","COL8A1"),
               rownames(seu))
DefaultAssay(seu) <- "SCT"
seu <- AddModuleScore(seu, features = list(x), assay = "SCT", name = "FAPModule")
DotPlot(seu, features = x, group.by = "MPsignatures_NN_res.0.5") + RotatedAxis()

avg <- AverageExpression(seu, assay = "SCT", group.by = "MPsignatures_NN_res.0.5", features = x)
pdf("./FAPModule_NMF_NN_res.0.5_2.pdf", useDingbats = FALSE, width = 2.5, height = 2.5)
pheatmap(avg$SCT, scale = "row", col = paletteContinuous("solarExtra"),
         cluster_rows = TRUE, fontsize_row = 5, fontsize_col = 5,
         cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

############################################################
# 8) MP overlap analyses (UpSet + FAP_module overlap bar)
############################################################
df <- read.csv("MP_genes.csv", stringsAsFactors = FALSE)

# Binary matrix of gene presence in MPs
binary_matrix <- df %>%
  distinct() %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from = MP, values_from = present, values_fill = list(present = 0)) %>%
  group_by(gene) %>% summarise(across(everything(), max), .groups = "drop")

pdf("./MP_UpSet.pdf", width = 6, height = 5.5)
upset(as.data.frame(binary_matrix %>% select(-gene)),
      nsets = ncol(binary_matrix) - 1, nintersects = 40,
      order.by = "freq", sets.bar.color = "steelblue", main.bar.color = "black")
dev.off()

FAP_module <- c("MGP","LUM","VCAN","F2R","OMD","LTBP2","FAP","FN1","COL1A2","CFH","PLXDC2",
                "THBS2","ITGBL1","PRELP","FBLN2","PDFRA","DCN","MMP2","PDGFD,","FXYD5",
                "POSTN","LAMA2","MEG3","COL3A1","TMSB10","CRTAC1","PCOLCE2","COL8A1")

df_MP5 <- subset.data.frame(df, df$MP == "MP5")
df_MP7 <- subset.data.frame(df, df$MP == "MP7")
intersect(df_MP5$gene, df_MP7$gene)

# Overlap counts per MP for FAP module
overlap_df <- binary_matrix %>% filter(gene %in% FAP_module)
binary_overlap <- overlap_df %>% select(-gene) %>% mutate(across(everything(), ~ ifelse(. > 0, 1, 0)))
overlap_counts <- binary_overlap %>% summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "MP", values_to = "OverlapCount") %>%
  arrange(desc(OverlapCount)) %>% filter(MP != "X")

pdf("./MP_FAPModuleOverlapp.pdf", width = 3, height = 3)
ggplot(overlap_counts, aes(x = reorder(MP, -OverlapCount), y = OverlapCount)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "FAP_module Gene Overlap with Each MP", x = "MP", y = "Overlapping Genes") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
