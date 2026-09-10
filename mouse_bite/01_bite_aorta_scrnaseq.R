############################################################
# (Control vs FapBiTE) — End-to-end analysis
# - Libraries deduplicated (kept only what’s used)
# - Original workflow preserved; added clear, concise comments
# - Notes:
#   * paletteDiscrete()/paletteContinuous() come from ArchR (loaded below).
#   * This script mixes loading fresh 10x data and later reading saved RDS.
#     That is intentional to mirror the original steps.
#   * Some object names (e.g., SR006237_FAP_SMC_dir as a Seurat object)
#     are unusual but left unchanged to preserve the original code flow.
############################################################

## ---- Libraries (unique & sufficient) ----
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(patchwork)
library(harmony)
library(ArchR)      # provides paletteDiscrete()/paletteContinuous()
library(ggsci)
library(ggpubr)     # used later for plotting styles
library(pheatmap)
library(readr)      # read_delim used for DE lists
library(clusterProfiler)  # for GO
library(org.Mm.eg.db)     # mouse OrgDb for GO
library(Augur)      # calculate_auc
options(future.globals.maxSize = 100000 * 1024^2) # ~100 GB

############################################################
# Preprocessing: Load 10x matrices, create objects, merge
############################################################

# Control (SMC-enriched)
SR006237_Control_SMC_dir  <- './Lavine_SR006237_10X/SR006237_Control_SMC/filtered_feature_bc_matrix/'
SR006237_Control_SMC.data <- Read10X(data.dir = SR006237_Control_SMC_dir)
SR006237_Control_SMC      <- CreateSeuratObject(counts = SR006237_Control_SMC.data)
SR006237_Control_SMC$condition <- "Control"

# FapBiTE (SMC-enriched)  -- object name preserved from original code
SR006237_FAP_SMC_dir      <- './Lavine_SR006237_10X/SR006237_FAP_SMC/filtered_feature_bc_matrix/'
SR006237_FAP_SMC_dir.data <- Read10X(data.dir = SR006237_FAP_SMC_dir)
SR006237_FAP_SMC_dir      <- CreateSeuratObject(counts = SR006237_FAP_SMC_dir.data)
SR006237_FAP_SMC_dir$condition <- "FapBiTE"

# Merge
sample <- merge(SR006237_Control_SMC, y = c(SR006237_FAP_SMC_dir))

############################################################
# QC: basic mitochondrial and feature/UMI filtering
############################################################

# Mito % and initial QC violin
sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^mt-")
VlnPlot(sample, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3, pt.size = 0, group.by = "condition")

# Filter cells (keep parameters from original)
sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)

# Post-filter QC violin by condition
VlnPlot(sample, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3, group.by = "condition")

# Load previously saved global object (as in original workflow)
sample <- readRDS("./global/v2/global.rds")

# Example panel: percent.mt across samples
pdf("./percent.mt.pdf", useDingbats = FALSE, width = 4, height = 3)
VlnPlot(sample, features = "percent.mt", ncol = 4, group.by = "orig.ident", pt.size = 0)
dev.off()

############################################################
# Normalization / integration-lite: SCTransform + PCA/UMAP/Neighbors/Clusters
############################################################
DefaultAssay(sample) <- 'RNA'
sample <- SCTransform(sample, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = TRUE)

sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 50, verbose = TRUE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:50)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:50)
sample <- FindClusters(sample, graph.name = "SCT_snn", algorithm = 3,
                       resolution = 0.3, verbose = FALSE)

# UMAP by clusters (split by treatment)
DimPlot(sample, reduction = 'umap', label = TRUE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.3",
        cols = paletteDiscrete(unique(sample$SCT_snn_res.0.3), set = "stallion"),
        split.by = "condition")

# Example feature (RNA)
FeaturePlot(sample, features = "Spp1")

# Cluster composition by condition
ggplot(sample@meta.data, aes(x = condition, fill = SCT_snn_res.0.3)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(sample$SCT_snn_res.0.3), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

############################################################
# DGE across clusters (SCT assay)
############################################################
Idents(sample) <- "SCT_snn_res.0.3"
DefaultAssay(sample) <- 'SCT'
rna.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./DE_SCT_snn_res.0.3.csv", quote = FALSE)

saveRDS(sample, "sample.rds")

# Load merged object for downstream annotation (per original)
sample <- readRDS("./global/merged.rds")
Idents(sample) <- "SCT_snn_res.0.3"

############################################################
# Cell-type labeling (string switch), remove “Junk”, recluster
############################################################
fun <- function(x) {
  if (x == "0") {"SMC"}
  else if (x == "1") {"Fibroblast"}
  else if (x == "2") {"Myeloid"}
  else if (x == "3") {"modSMC"}
  else if (x == "4") {"Junk"}
  else if (x == "5") {"Endothelium"}
  else if (x == "6") {"Pericyte"}
  else if (x == "7") {"Junk"}
  else if (x == "8") {"Myeloid"}
  else if (x == "9") {"TNKCell"}
  else if (x == "10") {"Fibroblast"}
  else if (x == "11") {"Junk"}
  else if (x == "12") {"Neutrophil"}
  else if (x == "13") {"Myeloid"}
  else if (x == "14") {"Junk"}
  else if (x == "15") {"Mesothelium"}
  else if (x == "16") {"Lymphatic"}
  else if (x == "17") {"Glia"}
  else if (x == "18") {"Proliferating"}
  else if (x == "19") {"Junk"}
}
sample$cell.type <- mapply(fun, sample$SCT_snn_res.0.3)

# Drop “Junk”, then re-run neighbors/clusters to refine
Idents(sample) <- "cell.type"
sample <- subset(sample, idents = "Junk", invert = TRUE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:50)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:50)
sample <- FindClusters(sample, graph.name = "SCT_snn",
                       algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

# UMAP & composition
DimPlot(sample, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "cell.type",
        cols = paletteDiscrete(unique(sample$cell.type), set = "stallion"))
ggplot(sample@meta.data, aes(x = condition, fill = cell.type)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(sample$cell.type), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

############################################################
# Per–cell type DGE: Control vs FapBiTE (loop)
############################################################
celltype <- unique(sample$cell.type)
for (cell in celltype) {
  Idents(sample) <- "cell.type"
  curr_subset <- subset(sample, idents = cell)
  Idents(curr_subset) <- "condition"
  DefaultAssay(curr_subset) <- 'SCT'
  rna.markers <- FindAllMarkers(curr_subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  write.csv(rna.markers,
            file = paste0("./global/v2/DE_lists/DE_IsovsFAPBiTE_", cell, ".csv"),
            quote = FALSE)
}

############################################################
# Summarize DE lists and visualize “significant” directions
############################################################
avg_log2FC_cutoff <- 0.25

# Read per-cell-type DGE tables, flag significance, and bind
Endothelium  <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Endothelium.csv", ",", show_col_types = FALSE)
Endothelium$cell <- "Endothelium"
Endothelium$sigpvalue <- ifelse(Endothelium$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Endothelium$sig <- ifelse(Endothelium$p_val_adj < 0.05 & Endothelium$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Fibroblast   <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Fibroblast.csv", ",", show_col_types = FALSE)
Fibroblast$cell <- "Fibroblast"
Fibroblast$sigpvalue <- ifelse(Fibroblast$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Fibroblast$sig <- ifelse(Fibroblast$p_val_adj < 0.05 & Fibroblast$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Glia <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Glia.csv", ",", show_col_types = FALSE)
Glia$cell <- "Glia"
Glia$sigpvalue <- ifelse(Glia$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Glia$sig <- ifelse(Glia$p_val_adj < 0.05 & Glia$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Lymphatic <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Lymphatic.csv", ",", show_col_types = FALSE)
Lymphatic$cell <- "Lymphatic"
Lymphatic$sigpvalue <- ifelse(Lymphatic$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Lymphatic$sig <- ifelse(Lymphatic$p_val_adj < 0.05 & Lymphatic$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Mesothelium <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Mesothelium.csv", ",", show_col_types = FALSE)
Mesothelium$cell <- "Mesothelium"
Mesothelium$sigpvalue <- ifelse(Mesothelium$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Mesothelium$sig <- ifelse(Mesothelium$p_val_adj < 0.05 & Mesothelium$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

modSMC <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_modSMC.csv", ",", show_col_types = FALSE)
modSMC$cell <- "modSMC"
modSMC$sigpvalue <- ifelse(modSMC$p_val_adj < 0.05, "p < 0.05","p > 0.05")
modSMC$sig <- ifelse(modSMC$p_val_adj < 0.05 & modSMC$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Myeloid <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Myeloid.csv", ",", show_col_types = FALSE)
Myeloid$cell <- "Myeloid"
Myeloid$sigpvalue <- ifelse(Myeloid$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Myeloid$sig <- ifelse(Myeloid$p_val_adj < 0.05 & Myeloid$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Neutrophil <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Neutrophil.csv", ",", show_col_types = FALSE)
Neutrophil$cell <- "Neutrophil"
Neutrophil$sigpvalue <- ifelse(Neutrophil$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Neutrophil$sig <- ifelse(Neutrophil$p_val_adj < 0.05 & Neutrophil$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Pericyte <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Pericyte.csv", ",", show_col_types = FALSE)
Pericyte$cell <- "Pericyte"
Pericyte$sigpvalue <- ifelse(Pericyte$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Pericyte$sig <- ifelse(Pericyte$p_val_adj < 0.05 & Pericyte$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

Proliferating <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_Proliferating.csv", ",", show_col_types = FALSE)
Proliferating$cell <- "Proliferating"
Proliferating$sigpvalue <- ifelse(Proliferating$p_val_adj < 0.05, "p < 0.05","p > 0.05")
Proliferating$sig <- ifelse(Proliferating$p_val_adj < 0.05 & Proliferating$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

SMC <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_SMC.csv", ",", show_col_types = FALSE)
SMC$cell <- "SMC"
SMC$sigpvalue <- ifelse(SMC$p_val_adj < 0.05, "p < 0.05","p > 0.05")
SMC$sig <- ifelse(SMC$p_val_adj < 0.05 & SMC$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

TNKCell <- read_delim("./global/v2/DE_lists/DE_IsovsFAPBiTE_TNKCell.csv", ",", show_col_types = FALSE)
TNKCell$cell <- "TNKCell"
TNKCell$sigpvalue <- ifelse(TNKCell$p_val_adj < 0.05, "p < 0.05","p > 0.05")
TNKCell$sig <- ifelse(TNKCell$p_val_adj < 0.05 & TNKCell$avg_log2FC > avg_log2FC_cutoff, "Significant","Not Significant")

# Bind all cell types; set factor orders
data <- data.frame(rbind(Endothelium,Fibroblast,Glia,Lymphatic,Mesothelium,
                         modSMC,Myeloid,Neutrophil,Pericyte,Proliferating,SMC,TNKCell))
data$cell <- factor(data$cell, levels = c("Endothelium","Fibroblast","Glia","Lymphatic",
                                          "Mesothelium","modSMC","Myeloid","Neutrophil",
                                          "Pericyte","Proliferating","SMC","TNKCell"))

# Flip sign so positive = Control-up, negative = BiTE-up (per original logic)
data$avg_log2FC <- ifelse(data$cluster == "Control", data$avg_log2FC, -1 * data$avg_log2FC)
table(data$cell, data$sig)

# Order cell types by count of significant genes
df_Count <- data %>% group_by(sig, cell) %>% dplyr::count()
df_Count <- data.frame(df_Count)
x <- df_Count[with(df_Count, order(n, decreasing = TRUE)), ] |>
     subset(sig == "Significant") |>
     (`[[`)("cell")
df_Count$cell <- factor(df_Count$cell, levels = x)
data$cell     <- factor(data$cell,     levels = x)

# Scatter of logFCs per cell type (colored by significance)
data %>%
  ggplot(aes(x = cell, y = avg_log2FC, fill = cell, color = sig)) +
  geom_jitter(size = 1, alpha = 0.5, position = position_jitter(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(limits = c(-1, 1)) +
  ggtitle("Pseudobulk DE WT_vs_KO") +
  xlab("Cell types") +
  scale_shape_manual(values = c(1, 1)) +
  scale_color_manual(values = c("grey", "red"))

# Re-derive directional label after sign flip
data$cluster <- ifelse(data$avg_log2FC > 0, "Control", "anti-FAP BiTE")

# Stacked bar of significant genes (Control-up vs BiTE-up) per cell type
pdf("./stack_DE.pdf", useDingbats = FALSE, width = 3.5, height = 2)
df_Count <- data %>% group_by(sig, cell) %>% dplyr::count(cluster)
df_Count <- data.frame(df_Count)
df_Count$cell <- factor(df_Count$cell)
df_Count <- dplyr::filter(df_Count, sig == "Significant")
ggplot(df_Count, aes(x = cell, y = n, fill = cluster)) +
  geom_col() +
  labs(x = "Cell Type", y = "#Genes", fill = "Direction") +
  scale_fill_manual(values = c("dodgerblue4","red3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()

############################################################
# AUGUR: perturbation separability per cell type
############################################################
sample_AUGR <- calculate_auc(sample, label_col = "condition", cell_type_col = "cell.type")
df_AUGR <- sample_AUGR$AUC
df_AUGR

############################################################
# Stromal focus: SMC/Pericyte/Fibroblast/modSMC
############################################################
Idents(sample) <- "cell.type"
stroma <- subset(sample, idents = c("SMC","Pericyte","Fibroblast","modSMC"))

stroma <- RunUMAP(stroma, reduction = "pca", dims = 1:50)
stroma <- FindNeighbors(stroma, reduction = "pca", dims = 1:50)
stroma <- FindClusters(stroma, graph.name = "SCT_snn", algorithm = 3,
                       resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

# Subcluster labels @ res=0.5 (later re-annotated as states)
DimPlot(stroma, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.5",
        cols = paletteDiscrete(unique(stroma$SCT_snn_res.0.5), set = "stallion"),
        ncol = 2)

# Example features
FeaturePlot(stroma, features = c("Myh11","Tnfrsf11b"))

# Composition by condition
ggplot(stroma@meta.data, aes(x = condition, fill = SCT_snn_res.0.5)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(stroma$SCT_snn_res.0.5), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# DGE across stroma subclusters
Idents(stroma) <- "SCT_snn_res.0.5"
DefaultAssay(stroma) <- 'SCT'
rna.markers <- FindAllMarkers(stroma, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./stroma/DE_SCT_snn_res.0.5.csv", quote = FALSE)

# State mapping (per original)
fun <- function(x) {
  if (x == "0") {"SMC1"}
  else if (x == "1") {"Fib1"}
  else if (x == "2") {"SMC2"}
  else if (x == "3") {"FMC2"}
  else if (x == "4") {"FMC1"}
  else if (x == "5") {"Pericyte"}
  else if (x == "6") {"Fib2"}
  else if (x == "7") {"Fib3"}
  else if (x == "8") {"Fib4"}
}
stroma$cell.state <- mapply(fun, stroma$SCT_snn_res.0.5)
stroma$cell.state <- factor(stroma$cell.state,
                            levels = c("Pericyte","SMC1","SMC2","FMC1","FMC2","Fib1","Fib2","Fib3","Fib4"))

# UMAP / composition by state
DimPlot(stroma, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "cell.state",
        cols = paletteDiscrete(unique(stroma$cell.state), set = "stallion"),
        ncol = 2)

ggplot(stroma@meta.data, aes(x = condition, fill = cell.state)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(stroma$cell.state), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Marker panel for states
pdf("./marker_dotplot.pdf", useDingbats = FALSE, width = 5.5, height = 3.5)
DotPlot(stroma, features = c("Rgs5","Kcnip4","Myh11","Acta2","Wnt16","Itih4",
                             "Tnfrsf11b","Fap","Abca8a","Fbln1")) + RotatedAxis()
dev.off()

# SMC module score (FMCz)
Idents(stroma) <- "cell.state"
SMCs <- subset(stroma, idents = c("Pericyte","SMC1","SMC2","FMC1","FMC2"))
DefaultAssay(SMCs) <- "SCT"
SMCs <- AddModuleScore(SMCs,
  features = list(c("Cnn1","Fap","Runx1","Tnfrsf11b","Dcn","Crip1","Spp1","Ccn2",
                    "Ltbp2","Col4a1","Ccn1","Atf3","Ccn5","Jun","Fos")),
  name = 'FMCz'
)
pdf("./human_FMCz.pdf", useDingbats = FALSE, width = 2.5, height = 3)
VlnPlot(SMCs, features = "FMCz1", group.by = "condition", pt.size = 0)
dev.off()

# FAP dotplot per state (re-load a saved object per original)
stroma <- readRDS("./stroma/stroma.rds")
pdf("./Fap_dotplot.pdf", useDingbats = FALSE, width = 3.5, height = 3.5)
DotPlot(stroma, features = "Fap", group.by = "cell.state") + RotatedAxis()
dev.off()

############################################################
# Cross-reference mapping to mouse reference (Quetermous/Cheng)
############################################################
mouse_paul <- readRDS("/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/analysis/final_analysis/Mouse_mapping/Quetermous_Mapping/Cheng_Zeb2/16wkhfdctl.rds")

stroma <- NormalizeData(stroma)  # prepare RNA assay for anchors

anchors <- FindTransferAnchors(reference = mouse_paul, query = stroma, dims = 1:30,
                               reference.reduction = "pca",
                               reference.assay = "RNA", query.assay = "RNA")

predictions <- TransferData(anchorset = anchors, refdata = mouse_paul$Genotype, dims = 1:30)
stroma <- AddMetaData(stroma, metadata = predictions)

stroma$predicted.id <- factor(stroma$predicted.id, levels = c("Ctl-tdt","Ctl-NOtdt"))

DimPlot(stroma, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "cell.state",
        cols = paletteDiscrete(unique(stroma$cell.state), set = "stallion"),
        split.by = "predicted.id")

saveRDS(stroma, "stroma_ref_mapped.rds")

############################################################
# tdT+ vs tdT- within stroma; per-condition DEG
############################################################
stroma <- readRDS("./stroma/stroma_ref_mapped.rds")

Idents(stroma) <- "predicted.id"
tdtPos <- subset(stroma, idents = "Ctl-tdt")
tdtNeg <- subset(stroma, idents = "Ctl-NOtdt")

Idents(tdtPos) <- "condition"
DefaultAssay(tdtPos) <- 'SCT'
rna.markers.tdtPos <- FindAllMarkers(tdtPos, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers.tdtPos, file = "./stroma/DE_tdTPos_condition.csv", quote = FALSE)

Idents(tdtNeg) <- "condition"
DefaultAssay(tdtNeg) <- 'SCT'
rna.markers.tdtNeg <- FindAllMarkers(tdtNeg, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers.tdtNeg, file = "./stroma/DE_tdtNeg_condition.csv", quote = FALSE)

# Filter significant & split by Control-up
rna.markers.tdtPos <- subset.data.frame(rna.markers.tdtPos, p_val_adj < 0.05)
rna.markers.tdtPos.Control <- subset.data.frame(rna.markers.tdtPos, cluster == "Control")

rna.markers.tdtNeg <- subset.data.frame(rna.markers.tdtNeg, p_val_adj < 0.05)
rna.markers.tdtNeg.Control <- subset.data.frame(rna.markers.tdtNeg, cluster == "Control")

############################################################
# Pathway enrichment (Mouse GO: Molecular Function)
############################################################
genes <- unique(rna.markers.tdtNeg.Control$gene)

mapping <- bitr(genes,
                fromType = "SYMBOL",
                toType   = "ENTREZID",
                OrgDb    = org.Mm.eg.db)
entrez <- mapping$ENTREZID

ego <- enrichGO(gene         = entrez,
                OrgDb        = org.Mm.eg.db,
                keyType      = "ENTREZID",
                ont          = "MF",        # MF/CC/BP
                pAdjustMethod= "BH",
                qvalueCutoff = 0.05)

pdf("./GO_MF_tdtNeg_downBiTE.pdf", useDingbats = FALSE, width = 5, height = 7)
dotplot(ego, showCategory = 20, font.size = 10)
dev.off()

############################################################
# Myeloid subset: recluster, annotate states, marker panels
############################################################
Idents(sample) <- "cell.type"
Myeloid <- subset(sample, idents = "Myeloid")

Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:50)
Myeloid <- FindNeighbors(Myeloid, reduction = "pca", dims = 1:50)
Myeloid <- FindClusters(Myeloid, graph.name = "SCT_snn", algorithm = 3,
                        resolution = c(0.1,0.2,0.3,0.4,0.5), verbose = FALSE)

DimPlot(Myeloid, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.5",
        cols = paletteDiscrete(unique(Myeloid$SCT_snn_res.0.5), set = "stallion"),
        ncol = 2)

ggplot(Myeloid@meta.data, aes(x = condition, fill = SCT_snn_res.0.5)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(Myeloid$SCT_snn_res.0.5), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# DGE per subcluster
Idents(Myeloid) <- "SCT_snn_res.0.5"
DefaultAssay(Myeloid) <- 'SCT'
rna.markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./Myeloid/DE_SCT_snn_res.0.5.csv", quote = FALSE)

# Remove cluster "3" (as in original), re-embed
Idents(Myeloid) <- "SCT_snn_res.0.5"
Myeloid <- subset(Myeloid, idents = "3", invert = TRUE)
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:50)

# Repeat plots & composition
DimPlot(Myeloid, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.5",
        cols = paletteDiscrete(unique(Myeloid$SCT_snn_res.0.5), set = "stallion"),
        ncol = 2)

ggplot(Myeloid@meta.data, aes(x = condition, fill = SCT_snn_res.0.5)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(Myeloid$SCT_snn_res.0.5), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# DGE again after filtering
Idents(Myeloid) <- "SCT_snn_res.0.5"
DefaultAssay(Myeloid) <- 'SCT'
rna.markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./myeloid/v2/DE_Myeloid_SCT_snn_res.0.5.csv", quote = FALSE)

# Map subclusters to states
fun <- function(x) {
  if (x == "0") {"Mac1"}
  else if (x == "1") {"Mac2"}
  else if (x == "2") {"Mac4"}
  else if (x == "4") {"Mac3"}
  else if (x == "5") {"cDC2"}
  else if (x == "6") {"Mono"}
  else if (x == "7") {"cDC1"}
  else if (x == "8") {"mDC"}
}
Myeloid$cell.state <- mapply(fun, Myeloid$SCT_snn_res.0.5)
Myeloid$cell.state <- factor(Myeloid$cell.state,
                             levels = c("Mono","Mac1","Mac2","Mac3","Mac4","cDC1","cDC2","mDC"))

# UMAP / composition by myeloid states
DimPlot(Myeloid, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "cell.state",
        cols = paletteDiscrete(unique(Myeloid$cell.state), set = "stallion"),
        ncol = 2)

ggplot(Myeloid@meta.data, aes(x = condition, fill = cell.state)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(Myeloid$cell.state), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# DGE across myeloid states
Idents(Myeloid) <- "cell.state"
DefaultAssay(Myeloid) <- 'SCT'
rna.markers <- FindAllMarkers(Myeloid, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./myeloid/v2/DE_Myeloid_cell.state.csv", quote = FALSE)

# Example feature across condition
FeaturePlot(Myeloid, features = "Gpnmb", split.by = "condition")

# Myeloid marker panel
pdf("./marker_myeloid_dotplot.pdf", useDingbats = FALSE, width = 6.3, height = 3.5)
DotPlot(Myeloid, features = c("Ly6c2","Plac8","Pparg","Apbb2","Fabp4","Fabp5",
                              "Gpnmb","Ctsl","Spp1","Cd163","F13a1","Clec9a",
                              "Cd209a","Ccr7")) + RotatedAxis()
dev.off()

# LAM-like program score
DefaultAssay(Myeloid) <- "SCT"
Myeloid <- AddModuleScore(Myeloid,
  features = list(c("Fabp4","Fabp5","Gpnmb","Ctsl","Spp1")),
  name = 'LAMz'
)

# Quick dot and violin panels (as in original)
pdf("./LAMz1_myeloid_dotplot.pdf", useDingbats = FALSE, width = 3.3, height = 3.5)
DotPlot(Myeloid, features = "LAMz1") + RotatedAxis()
dev.off()

# Example violin for a single gene by condition
VlnPlot(Myeloid, features = "Ctsl", group.by = "condition", pt.size = 0)
