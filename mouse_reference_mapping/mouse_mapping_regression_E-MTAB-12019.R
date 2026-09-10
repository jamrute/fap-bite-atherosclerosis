#########################################
# Libraries
#########################################
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ArchR)
library(ggsci)

###### Preprocessing Steps ######

#########################################
# Load 10x datasets and create Seurat objects
# Each object gets a 'condition' label for downstream integration/batch correction
#########################################

# ApoBInjected (baseline)
ApoBInjected_dir <- './data/ApoBInjected/'
ApoBInjected.data <- Read10X(data.dir = ApoBInjected_dir)              # reads 10x matrix
ApoBInjected <- CreateSeuratObject(counts = ApoBInjected.data)         # Seurat object
ApoBInjected$condition <- "ApoBInjected"                               # annotate condition

# ApoBInjected at 12 weeks
ApoBInjected_12W_dir <- './data/ApoBInjected_12W/'
ApoBInjected_12W.data <- Read10X(data.dir = ApoBInjected_12W_dir)
ApoBInjected_12W <- CreateSeuratObject(counts = ApoBInjected_12W.data)
ApoBInjected_12W$condition <- "ApoBInjected_12W"

# High-fat diet control (baseline)
HFDControl_dir <- './data/HFDControl/'
HFDControl.data <- Read10X(data.dir = HFDControl_dir)
HFDControl <- CreateSeuratObject(counts = HFDControl.data)
HFDControl$condition <- "HFDControl"

# High-fat diet control at 12 weeks
HFDControl_12W_dir <- './data/HFDControl_12W/'
HFDControl_12W.data <- Read10X(data.dir = HFDControl_12W_dir)
HFDControl_12W <- CreateSeuratObject(counts = HFDControl_12W.data)
HFDControl_12W$condition <- "HFDControl_12W"

# Merge all samples into a single Seurat object
# (cell barcodes are automatically made unique internally)
sample <- merge(ApoBInjected, y = c(ApoBInjected_12W, HFDControl, HFDControl_12W))

#########################################
# QC features and initial visualization
#########################################

# Compute mitochondrial percentage
# NOTE: for mouse use "^mt-", for human use "^MT-"
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")

# Visualize basic QC metrics by condition (no points for speed/clarity)
VlnPlot(sample,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0, group.by = "condition")

# Basic QC filtering thresholds (tune as needed)
sample <- subset(sample,
                 subset = nFeature_RNA > 200 &
                          nFeature_RNA < 6000 &
                          percent.mt   < 10 &
                          nCount_RNA   < 35000)

# Re-plot QC distributions after filtering
VlnPlot(sample,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, group.by = "condition")

#########################################
# Normalization (SCTransform), dimensionality reduction, and clustering
#########################################

# Work in RNA assay prior to SCT
DefaultAssay(sample) <- 'RNA'

# SCTransform with regression of mt% and total RNA counts (depth)
sample <- SCTransform(sample,
                      vars.to.regress = c("percent.mt", "nCount_RNA"),
                      verbose = TRUE)

# PCA on variable features; compute 100 PCs (you will likely use fewer downstream)
sample <- RunPCA(sample,
                 features = VariableFeatures(object = sample),
                 npcs = 100, verbose = TRUE)

# Harmony batch correction using 'condition'
# - reduction.save = "harmony" stores Harmony embedding
# - assay.use = "SCT" is appropriate post-SCTransform
sample <- RunHarmony(sample, c("condition"),
                     reduction = "pca",
                     reduction.save = "harmony",
                     assay.use = "SCT")

# UMAP/Neighbors/Clustering on Harmony embedding (dims 1:50 used here)
sample <- RunUMAP(sample, reduction = "harmony", dims = 1:50)
sample <- FindNeighbors(sample, reduction = "harmony", dims = 1:50)

# Leiden clustering (algorithm = 3). Multiple resolutions recorded simultaneously.
sample <- FindClusters(sample,
                       graph.name = "SCT_snn",
                       algorithm = 3,
                       resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                       verbose = FALSE)

# UMAP colored by a specific resolution (0.2). Uses a palette from your env (paletteDiscrete).
# If paletteDiscrete is not defined, consider replacing with a ggplot2/ggsci palette.
DimPlot(sample, reduction = 'umap', label = TRUE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.2",
        cols = paletteDiscrete(unique(sample$SCT_snn_res.0.2), set = "stallion"))

#########################################
# Differential Gene Expression (cluster markers on SCT assay)
#########################################

# Set identities to the chosen clustering level
Idents(sample) <- "SCT_snn_res.0.2"
DefaultAssay(sample) <- 'SCT'

# Positive markers per cluster with minimal prevalence and effect size filters
rna.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

# Export DGE table
write.csv(rna.markers, file = "./DE_SCT_snn_res.0.2.csv", quote = FALSE)

# Save the whole object (post-clustering)
saveRDS(sample, "global.rds")

#########################################
# Biological subsetting (edit cluster IDs to match your biology)
#########################################

# Re-affirm identities and order the condition factor for consistent plotting
Idents(sample) <- "SCT_snn_res.0.2"
sample$condition <- factor(sample$condition,
                           levels = c("HFDControl","HFDControl_12W","ApoBInjected","ApoBInjected_12W"))

# Example compartment splits (verify these cluster IDs!)
endothelium <- subset(sample, idents = c("6"))
stroma      <- subset(sample, idents = c("0","1","2","3","5","7","8","9"))
myeloid     <- subset(sample, idents = c("4"))

# Example marker visualization in stroma by condition
DotPlot(stroma, features = c("Fap","Runx1","Sox9","Myh11"), group.by = "condition") + RotatedAxis()
DotPlot(stroma, features = c("Spp1","Lgals3"), group.by = "condition") + RotatedAxis()

# Quick check of condition composition in myeloid subset
unique(myeloid$condition)

#########################################
# Condition-restricted DE: myeloid (SCT)
#########################################

Idents(myeloid) <- "condition"
myeloid2 <- subset(myeloid, idents = c("HFDControl_12W","ApoBInjected_12W"))
DefaultAssay(myeloid2) <- 'SCT'

rna.markers <- FindAllMarkers(myeloid2, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(rna.markers, file = "./DE_myeloid_ApoBInjected_12W_vs_HFDControl_12W.csv", quote = FALSE)

#########################################
# Condition-restricted DE: stroma (SCT)
#########################################

Idents(stroma) <- "condition"
stroma2 <- subset(stroma, idents = c("HFDControl_12W","ApoBInjected_12W"))
DefaultAssay(stroma2) <- 'SCT'

rna.markers <- FindAllMarkers(stroma2, min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(rna.markers, file = "./DE_stroma_ApoBInjected_12W_vs_HFDControl_12W.csv", quote = FALSE)

# (Optional) re-save the global object again after subsets/DE (duplicate but kept as requested)
saveRDS(sample, "global.rds")

#########################################
# Additional libraries for downstream analyses/plots
#########################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(dplyr)
library(sctransform)
library(pheatmap)
library(Matrix)
library(RColorBrewer)
library(scales)
library(data.table)
library(stats)
library(ggsci)
library(ArchR)
library(biomaRt)
library(nichenetr)
library(tidyverse)

#########################################
# Load annotated human myeloid reference (Seurat object)
#########################################
ref <- readRDS("./myeloid_annotated.rds")

#########################################
# Mouse → Human gene symbol conversion and Seurat object creation
# Assumes you have a function convert_mouse_to_human_symbols() in your environment.
#########################################
mydata <- myeloid
mouse_rna_matrix <- mydata@assays[["RNA"]]@counts

# Convert gene names to human symbols; drop rows/cols with NAs post-conversion
rownames(mouse_rna_matrix) <- mouse_rna_matrix %>% rownames() %>% convert_mouse_to_human_symbols()
mouse_rna_matrix <- mouse_rna_matrix %>% .[!is.na(rownames(mouse_rna_matrix)), !is.na(colnames(mouse_rna_matrix))]

# Build a new Seurat object with converted counts and original metadata
mouse_coronary_new <- CreateSeuratObject(counts = mouse_rna_matrix)
mouse_coronary_new <- AddMetaData(object = mouse_coronary_new, metadata = mydata@meta.data)

# Keep only genes present in the reference for mapping
DefaultAssay(mouse_coronary_new) <- "RNA"
DefaultAssay(ref) <- "RNA"
mouse_coronary_new <- subset(mouse_coronary_new, features = rownames(ref))

#########################################
# Process the mouse (now human symbols) dataset prior to SCT
#########################################
DefaultAssay(mouse_coronary_new) <- "RNA"
mouse_coronary_new <- NormalizeData(mouse_coronary_new)
mouse_coronary_new <- FindVariableFeatures(mouse_coronary_new, selection.method = "vst", nfeatures = 3000)
mouse_coronary_new <- ScaleData(mouse_coronary_new)

# SCTransform for improved normalization/integration
mouse_coronary_new <- SCTransform(mouse_coronary_new, verbose = FALSE)

# Save intermediate
saveRDS(mouse_coronary_new, "mouse_to_human_normalized.rds")

#########################################
# Reference mapping
# - RunSPCA presumes an existing 'RNA_snn' graph in the reference (as in Seurat v5 workflows)
# - Ensure your 'ref' has the reductions/graphs requested below (spca, rna.umap)
#########################################
human_coronary <- ref
human_coronary <- RunSPCA(human_coronary, graph = "RNA_snn")

DefaultAssay(human_coronary) <- "SCT"
DefaultAssay(mouse_coronary_new) <- "SCT"

anchors <- FindTransferAnchors(
  reference = human_coronary,
  query = mouse_coronary_new,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  reference.assay = "SCT"
)

mouse_coronary_new <- MapQuery(
  anchorset = anchors,
  query = mouse_coronary_new,
  reference = human_coronary,
  refdata = list(
    celltype = "cell.state",  # metadata column in reference with labels
    predicted_ADT = "ADT"     # optional: predicted ADT if present in reference
  ),
  reference.reduction = "spca",
  reduction.model = "rna.umap" # must exist in reference
)

# UMAP plot of mapped query labeled by predicted cell type
DimPlot(mouse_coronary_new, reduction = "ref.umap", group.by = "predicted.celltype", label.size = 4,
        cols = paletteDiscrete(unique(human_coronary$cell.state), set = "stallion"),
        label = FALSE, split.by = "condition")

# Stacked bar chart of predicted cell types per condition
ggplot(mouse_coronary_new@meta.data, aes(x = condition, fill = predicted.celltype)) +
  geom_bar(position = "fill") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = paletteDiscrete(unique(human_coronary$cell.state), set = "stallion")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#########################################
# DE across conditions on the mapped object (RNA assay)
#########################################
Idents(mouse_coronary_new) <- "condition"
DefaultAssay(mouse_coronary_new) <- 'RNA'
rna.markers <- FindAllMarkers(mouse_coronary_new, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./DE_condition.csv", quote = FALSE)

# Split DEG table into "progression" vs "regression" clusters
prog <- rna.markers %>% filter(rna.markers$cluster == "progression")
reg  <- rna.markers %>% filter(rna.markers$cluster == "regression")

# Select top-200 genes by avg_log2FC for each program
prog %>%
  group_by(cluster) %>%
  top_n(n = 200, wt = avg_log2FC) -> prog_top10

reg %>%
  group_by(cluster) %>%
  top_n(n = 200, wt = avg_log2FC) -> reg_top10

#########################################
# Project progression/regression signatures into the reference
# (Manual z-scoring per gene then averaging; kept as in your code)
#########################################
DefaultAssay(ref) <- "RNA"
expdata <- GetAssayData(ref)
Pop1 <- prog_top10$gene
pops <- list(Pop1)

# Progression z-scores (manual approach preserved)
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
ref@meta.data$progression <- z_scores[1,]
FeaturePlot(object = ref, features = "progression", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0, 0.5))

# Regression z-scores
DefaultAssay(ref) <- "RNA"
expdata <- GetAssayData(ref)
Pop1 <- reg_top10$gene
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
ref@meta.data$regression <- z_scores[1,]
FeaturePlot(object = ref, features = "regression", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0, 1))

# Ratio and difference between programs (offsets guard against negatives/zeros)
ref$reg_to_prog <- (ref$regression + (min(ref$regression) * -1)) /
                   (ref$progression + (min(ref$progression) * -1))
ref$reg_minus_prog <- ref$regression - ref$progression

# Distributions by cell state
VlnPlot(ref, features = "reg_to_prog",   group.by = "cell.state", sort = TRUE, pt.size = 0, log = TRUE)
VlnPlot(ref, features = "reg_minus_prog", group.by = "cell.state", sort = TRUE, pt.size = 0, log = TRUE)

#########################################
# Foam/Niche signatures (manual z-scores preserved)
#########################################

# Mouse myeloid FoamNiche on SCT assay
DefaultAssay(myeloid) <- "SCT"
expdata <- GetAssayData(myeloid)
Pop1 <- c("Spp1","Cd74","Ftl1","Ctsd","Ctsb","H2-Aa","Apoe","Timp1","Laptm5","Fth1","C1qc","Psap")
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
myeloid@meta.data$FoamNiche <- z_scores[1,]

# Plot FoamNiche by condition (PDF output)
pdf("./FoamNiche_mouse_regression.pdf", useDingbats = FALSE, width = 4.3, height = 2.8)
DotPlot(myeloid, features = "FoamNiche", group.by = "condition") + RotatedAxis()
dev.off()

# Niche 5 signature in CITE-seq data (human symbols; SCT assay)
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

# End of script
