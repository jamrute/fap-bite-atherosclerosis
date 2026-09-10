############################################################
# Myeloid-focused analysis & foam-cell mapping (CITE-seq)
# - Libraries deduplicated (kept only what’s used)
# - Original code flow preserved; added concise comments
# - paletteDiscrete/paletteContinuous come from ArchR
############################################################

## ---- Libraries (unique & sufficient) ----
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(patchwork)
library(harmony)
library(ArchR)     # for paletteDiscrete/paletteContinuous
library(ggsci)
library(readr)
library(pheatmap)
library(progeny)
library(tidyr)
library(tibble)
library(Nebulosa)
library(viridis)

############################################################
# Load global object and subset to Myeloid
############################################################
global <- readRDS("./integrated_annotated.rds")

Idents(global) <- "cell.type"
levels(global)

sample <- subset(global, idents = c("Myeloid"))

############################################################
# RNA workflow: Normalize/Scale/Variable features/PCA/Harmony/UMAP/Neighbors/Clusters
############################################################
DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample)
all.genes <- rownames(sample)
sample <- ScaleData(sample, features = all.genes)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 100, verbose = TRUE)
sample <- RunHarmony(sample, c("sampleID"), reduction = "pca", reduction.save = "harmony_rna", assay.use = "RNA")
sample <- RunUMAP(sample, reduction = "harmony_rna", dims = 1:50,
                  reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', assay = "RNA")
sample <- FindNeighbors(sample, reduction = "harmony_rna", dims = 1:50)
sample <- FindClusters(sample, graph.name = "RNA_snn", algorithm = 3,
                       resolution = c(0.1, 0.2, 0.3), verbose = FALSE)

# ADT normalization (CLR) & scaling for later protein-level analyses
DefaultAssay(sample) <- "ADT"
sample <- NormalizeData(sample, normalization.method = 'CLR', margin = 2) %>% ScaleData()

# UMAP by RNA clusters
DimPlot(sample, reduction = 'rna.umap', group.by = 'RNA_snn_res.0.3',
        label.size = 4, label = TRUE,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.3), set = "stallion"))

############################################################
# DEG across RNA_snn_res.0.3 (RNA)
############################################################
DefaultAssay(sample) <- 'RNA'
Idents(sample) <- "RNA_snn_res.0.3"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./v1/v1_DE_RNA_snn_res.0.3.csv", quote = FALSE)

# Keep main clusters 0–9 only and re-embed/recluster at multiple resolutions
Idents(sample) <- "RNA_snn_res.0.3"
levels(sample)
sample <- subset(sample, idents = c("0","1","2","3","4","5","6","7","8","9"))

sample <- RunUMAP(sample, reduction = "harmony_rna", dims = 1:50,
                  reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', assay = "RNA", return.model = TRUE)
sample <- FindNeighbors(sample, reduction = "harmony_rna", dims = 1:50)
sample <- FindClusters(sample, graph.name = "RNA_snn", algorithm = 3,
                       resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

DimPlot(sample, reduction = 'rna.umap', group.by = 'RNA_snn_res.0.5',
        label.size = 4, label = TRUE,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.5), set = "stallion"))

# DEG across RNA_snn_res.0.5 (RNA)
DefaultAssay(sample) <- 'RNA'
Idents(sample) <- "RNA_snn_res.0.5"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./v2/v2_DE_RNA_snn_res.0.5.csv", quote = FALSE)

# DEG across RNA_snn_res.0.3 (ADT)
DefaultAssay(sample) <- 'ADT'
Idents(sample) <- "RNA_snn_res.0.3"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./v2_DE_ADT_snn_res.0.3.csv", quote = FALSE)

saveRDS(sample, "./v2/v2.rds")

# Remove cluster "10" if present and re-embed/recluster
Idents(sample) <- "RNA_snn_res.0.5"
sample <- subset(sample, idents = "10", invert = TRUE)
sample <- RunUMAP(sample, reduction = "harmony_rna", dims = 1:50,
                  reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', assay = "RNA", return.model = TRUE)
sample <- FindNeighbors(sample, reduction = "harmony_rna", dims = 1:50)
sample <- FindClusters(sample, graph.name = "RNA_snn", algorithm = 3,
                       resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

DimPlot(sample, reduction = 'rna.umap', group.by = 'RNA_snn_res.0.5',
        label.size = 4, label = TRUE,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.5), set = "stallion"))

# Final DEG after refinement (RNA)
DefaultAssay(sample) <- 'RNA'
Idents(sample) <- "RNA_snn_res.0.5"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./v3_DE_RNA_snn_res.0.5.csv", quote = FALSE)

############################################################
# Annotate RNA_snn_res.0.5 to biologically meaningful states
############################################################
fun <- function(x) {
  if (x == "0") {"Mac1"} 
  else if (x == "1") {"Mac2"}
  else if (x == "2") {"Mac6"}
  else if (x == "3") {"cMono"}
  else if (x == "4") {"cDC2"}
  else if (x == "5") {"Mac3"}
  else if (x == "6") {"Mac4"}
  else if (x == "7") {"Mac7"}
  else if (x == "8") {"Mac5"}
  else if (x == "9") {"ncMono"}
  else if (x == "10") {"Prolif"}
  else if (x == "11") {"cDC1"}
}
sample$cell.state <- mapply(fun, sample$RNA_snn_res.0.5)
sample$cell.state <- factor(sample$cell.state,
                            levels = c("cMono","ncMono","Mac1","Mac2","Mac3","Mac4","Mac5","Mac6","Mac7","cDC1","cDC2","Prolif"))

DimPlot(sample, reduction = 'rna.umap', group.by = 'cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(sample$cell.state), set = "stallion"))

# DEG by annotated state (RNA & ADT)
DefaultAssay(sample) <- 'RNA'
Idents(sample) <- "cell.state"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./DE_RNA_cell.state.csv", quote = FALSE)

DefaultAssay(sample) <- 'ADT'
Idents(sample) <- "cell.state"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file = "./DE_ADT_cell.state.csv", quote = FALSE)

saveRDS(sample, "./myeloid_annotated.rds")

############################################################
# Custom FOAM score (RNA z-score module) & protein features
############################################################
DefaultAssay(sample) <- "RNA"
expdata <- GetAssayData(sample)
Pop1 <- c("FABP4","FABP5","GPNMB","CTSL","HTRA1","PLD3")
pops <- list(Pop1)
z_scores <- NULL
for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores, (av + colSums(geneExp) / length(zz)))
}
sample@meta.data$FOAM <- z_scores[1,]
FeaturePlot(sample, features = "FOAM", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,2))

DefaultAssay(sample) <- "ADT"
FeaturePlot(sample, reduction = 'rna.umap', features = "CD86.1") +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,2))

############################################################
# Niche 5 signature score on SCT assay
############################################################
DefaultAssay(sample) <- "SCT"
expdata <- GetAssayData(sample)
Pop1 <- c("SPP1","CD74","FTL","CTSD","CTSB","HLA-DRA","APOE","TIMP1","LAPTM5","FTH1","C1QC","PSAP")
pops <- list(Pop1)
z_scores <- NULL
for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores, (av + colSums(geneExp) / length(zz)))
}
sample@meta.data$Niche5 <- z_scores[1,]

FeaturePlot(sample, features = "Niche5", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,1))

# Example GEX/ADT features
DefaultAssay(sample) <- "SCT"
FeaturePlot(sample, features = "LPL", reduction = 'rna.umap') +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,2))

DefaultAssay(sample) <- "ADT"
FeaturePlot(sample, features = "CCR2.1", reduction = 'rna.umap') +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,2))

############################################################
# PROGENy pathway activity (per cell.state)
############################################################
Idents(sample) <- "cell.state"
CellsClusters <- data.frame(Cell = names(Idents(sample)),
                            CellType = as.character(Idents(sample)),
                            stringsAsFactors = FALSE)

# Compute PROGENy scores (new assay) and scale them
sample <- progeny(sample, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
sample <- Seurat::ScaleData(sample, assay = "progeny")

# Summarize pathway activity by cell.state and heatmap
progeny_scores_df <- as.data.frame(t(GetAssayData(sample, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) %>%
  inner_join(CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>% spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue","white","red"))(paletteLength)
progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, length.out = ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength,
                       max(summarized_progeny_scores_df), length.out = floor(paletteLength/2)))
progeny_hmap <- pheatmap(t(summarized_progeny_scores_df[,-1]),
                         fontsize = 14, fontsize_row = 10,
                         color = myColor, breaks = progenyBreaks,
                         main = "PROGENy (500)", angle_col = 45,
                         border_color = NA, scale = "row")

############################################################
# Load finalized objects and visualize mapping results
############################################################
sample <- readRDS("./myeloid_annotated.rds")
foam   <- readRDS("./Foam_Cell_10x.rds")

DimPlot(sample, reduction = 'rna.umap', group.by = 'cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(sample$cell.state), set = "stallion"))

DimPlot(foam, reduction = 'ref.umap', group.by = 'predicted.cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(sample$cell.state), set = "stallion"))

foam$data <- "InVivo"
pdf("./foam_mapping_ridge.pdf", useDingbats = FALSE, width = 4, height = 2)
RidgePlot(foam, features = "predicted.cell.state.score", group.by = "data")
dev.off()

ggplot(foam@meta.data, aes(x = treatment, fill = predicted.cell.state)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = paletteDiscrete(unique(sample$cell.state), set = "stallion")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

############################################################
# Foam niche score on foam object; quick plot
############################################################
DefaultAssay(foam) <- "RNA"
expdata <- GetAssayData(foam)
Pop1 <- c("SPP1","CD74","FTL","CTSD","CTSB","HLA-DRA","APOE","TIMP1","LAPTM5","FTH1","C1QC","PSAP")
pops <- list(Pop1)
z_scores <- NULL
for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores, (av + colSums(geneExp) / length(zz)))
}
foam@meta.data$Foamz <- z_scores[1,]

DotPlot(foam, features = "Foamz", group.by = "treatment") + RotatedAxis()

pdf("./FoamNiche_score_invitro.pdf", useDingbats = FALSE, width = 4.3, height = 2.8)
DotPlot(myeloid, features = "Foamz", group.by = "treatment") + RotatedAxis()
dev.off()

# Protein features and density
DefaultAssay(sample) <- "ADT"
FeaturePlot(sample, reduction = 'rna.umap', features = "TFRC.1") +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,2))
plot_density(sample, features = "CD276.1", reduction = "rna.umap")

# FOAM score on sample (repeat for convenience)
DefaultAssay(sample) <- "RNA"
expdata <- GetAssayData(sample)
Pop1 <- c("FABP4","FABP5","GPNMB","CTSL","HTRA1","PLD3")
pops <- list(Pop1)
z_scores <- NULL
for (i in 1:length(pops)) {
  genes <- pops[[i]]
  zz <- which(tolower(rownames(expdata)) %in% tolower(genes))
  av <- numeric(ncol(expdata))
  geneExp <- as.matrix(expdata[zz, ])
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  z_scores <- rbind(z_scores, (av + colSums(geneExp) / length(zz)))
}
sample@meta.data$FOAM <- z_scores[1,]
FeaturePlot(sample, features = "FOAM", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,2))

DotPlot(sample, features = c("FABP4","FABP5","GPNMB","CTSL","SPP1"), group.by = "cell.state") + RotatedAxis()

pdf("./foam_genes.pdf", useDingbats = FALSE, width = 4.3, height = 3.5)
DotPlot(sample, features = c("FABP4","FABP5","GPNMB","CTSL","SPP1"), group.by = "cell.state") + RotatedAxis()
dev.off()

############################################################
# GWAS overlap pipeline (using precomputed DE)
############################################################
sample <- readRDS("./final/myeloid_annotated.rds")

d <- read.csv("./final/myeloid_DE_RNA_cell.state.csv")
d <- d %>% filter(avg_log2FC > 0.58, p_val_adj < 0.05)

d %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
top10

# GEX heatmap of top markers per state
DefaultAssay(sample) <- "SCT"
pdf("./GEX_marker_dotplot.pdf", useDingbats = FALSE, width = 16, height = 4)
DotPlot(sample, features = unique(top10$gene), group.by = "cell.state",
        col.min = 0, cols = c("lightgrey","red")) + RotatedAxis()
dev.off()

# Example ADT features
DefaultAssay(sample) <- "ADT"
FeaturePlot(sample, reduction = 'rna.umap', features = "LILRA4.1") +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,2))

# Overlap DGE with CAD GWAS list
DGE <- d
DGE <- filter(DGE, p_val_adj < 0.05)
diff_exp_genes <- unique(DGE$gene)
length(diff_exp_genes)

GWAS <- read_delim("./cad_gwas_genes.txt", ",", escape_double = FALSE, trim_ws = TRUE)
gwas_genes <- list(GWAS$gene)[[1]]
gwas_genes

intersecting_genes <- Reduce(intersect, list(diff_exp_genes, gwas_genes))
length(intersecting_genes)

DGE_subset <- filter(DGE, gene %in% intersecting_genes)
DGE_subset

cad_gwas_overlap <- data.frame(sort(table(DGE_subset$cluster), decreasing = TRUE))
cad_gwas_overlap

# Heatmap of GWAS-overlap genes by cell.state (single-cell averages)
avg <- AverageExpression(sample, assay = "SCT", group.by = "cell.state", features = DGE_subset$gene)
avg

pdf("./MyeloidGWASGenes_celstate_heatmap.pdf", useDingbats = FALSE, width = 2.5, height = 2)
pheatmap(avg$SCT, scale = "row", col = paletteContinuous("solarExtra"),
         cluster_rows = TRUE, fontsize_row = 5, fontsize_col = 5,
         cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

############################################################
# Visium spatial: feature & module scores
############################################################
visium <- readRDS("/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/Spatial/Visium_FFPE/analysis/seurat/integrated.rds")

SpatialFeaturePlot(object = visium, features = "LIPA", images = "Visium_T1077") +
  scale_fill_gradientn(colors = viridis::inferno(256))

avg <- AverageExpression(visium, assay = "SCT", group.by = "MPsignatures_NN_res.0.5", features = DGE_subset$gene)
avg

pdf("./MyeloidGWASGenes_Space_heatmap.pdf", useDingbats = FALSE, width = 2.5, height = 2)
pheatmap(avg$SCT, scale = "row", col = paletteContinuous("solarExtra"),
         cluster_rows = TRUE, fontsize_row = 5, fontsize_col = 5,
         cluster_cols = TRUE, legend = TRUE, border_color = NA)
dev.off()

# Mac4 gene-set module on Visium
d <- read.csv("./final/myeloid_DE_RNA_cell.state.csv") %>%
     filter(avg_log2FC > 0.58, p_val_adj < 0.05, cluster == 'Mac4')
d %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(visium) <- "SCT"
visium <- AddModuleScore(visium, features = list(top10$gene), assay = "SCT", name = "Mac4GS")

SpatialFeaturePlot(object = visium, features = "FAPModule1", images = "Visium_T1077") +
  scale_fill_gradientn(colors = viridis::inferno(256))
SpatialFeaturePlot(object = visium, features = "Mac4GS1", images = "Visium_T1077") +
  scale_fill_gradientn(colors = viridis::inferno(256))

############################################################
# In vitro foam (foam) clustering & scores
############################################################
DimPlot(foam, reduction = 'umap.harmony', group.by = 'RNA_snn_res.0.2',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(foam$RNA_snn_res.0.2), set = "stallion"))

ggplot(foam@meta.data, aes(x = treatment, fill = RNA_snn_res.0.2)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = paletteDiscrete(unique(foam$RNA_snn_res.0.2), set = "stallion")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

DefaultAssay(foam) <- "RNA"
foam <- AddModuleScore(foam, features = list(top10$gene), assay = "RNA", name = "Mac4GS")

plot_density(foam, features = "Mac4GS1", reduction = "umap.harmony")
top10

# Annotate foam RNA_snn_res.0.2 to IVMac* states
fun <- function(x) {
  if (x == "1") {"IVMac1"} 
  else if (x == "2") {"IVMac2"}
  else if (x == "3") {"IVMac3"}
  else if (x == "4") {"IVMac4"}
  else if (x == "5") {"IVMac5"}
  else if (x == "6") {"IVMac6"}
  else if (x == "7") {"IVMac7"}
}
foam$cell.state <- mapply(fun, foam$RNA_snn_res.0.2)

DimPlot(foam, reduction = 'umap.harmony', group.by = 'cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(foam$cell.state), set = "stallion"))

ggplot(foam@meta.data, aes(x = treatment, fill = cell.state)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = paletteDiscrete(unique(foam$cell.state), set = "stallion")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

pdf("./IVMac_cellstate_marker_dotplot.pdf", useDingbats = FALSE, width = 7, height = 3)
DotPlot(foam, features = c("FABP4","GPNMB","NPTX1","TYMS","MYBL2","FBXO5","STARD4",
                           "MSMO1","DERL3","PTTG1","PRTN3","CDC20","TOP2A",
                           "CDK1","BIRC5","MX1","ISG15","OAS2","MIF","H2AFZ","STMN1"),
        group.by = "cell.state") + RotatedAxis()
dev.off()

############################################################
# Selected features and LAM-like score in myeloid
############################################################
sample <- readRDS("./final/myeloid_annotated.rds")

DefaultAssay(sample) <- "SCT"
FeaturePlot(sample, features = "SIRPA")
DotPlot(sample, features = "SIRPA", group.by = "cell.state")

LAM_gene_signature <- c("CD63","CD68","CSTB","FABP5","GPNMB","LGALS3","LIPA","PLD3","PLIN2",
                        "SPP1","TREM2","CTSD","LPL")
DefaultAssay(sample) <- "SCT"
sample <- AddModuleScore(sample, features = list(LAM_gene_signature), assay = "SCT", name = "LAM")

pdf("./LAM_score_myeloid.pdf", useDingbats = FALSE, width = 3, height = 3)
RidgePlot(sample, features = "LAM1", group.by = "cell.state",
          cols = paletteDiscrete(unique(sample$cell.state), set = "stallion")) + NoLegend()
dev.off()
