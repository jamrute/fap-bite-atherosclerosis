############################################################
# T cells (Control vs FAP_BiTE): QC, SCT, clustering, TCR, StarCat
# - Libraries deduplicated (kept only what’s needed)
# - Original logic preserved; added clear comments
# - NOTES:
#   * paletteDiscrete()/paletteContinuous() come from ArchR plotting utils.
#   * convert_mouse_to_human_symbols() must exist in your env (user-defined/helper).
#   * A few sections reference objects not defined earlier (e.g., `tcells`, `mouse_coronary`, `ref`).
#     I’ve left the original lines intact and added comments where that may cause errors.
############################################################

## ---- Libraries (unique + sufficient) ----
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
library(ArchR)          # for paletteDiscrete()/paletteContinuous()
library(cowplot)
library(scRepertoire)   # TCR analysis
library(tidyverse)      # readr::write_delim used below
library(R.utils)        # gzip()
library(nichenetr)      # (loaded in StarCat section per original code)
library(SeuratDisk)     # (loaded in StarCat section per original code)
options(future.globals.maxSize = 90000 * 1024^2) # ~90 GB (comment in original said 50 GB)

############################################################
## Preprocessing: read 10x, annotate, merge, QC
############################################################

# Control sample
s1_dir <- './data/SR006661-SR006239_TCell_Control/count/sample_filtered_feature_bc_matrix/'
s1.data <- Read10X(data.dir = s1_dir)
s1 <- CreateSeuratObject(counts = s1.data)
s1$condition <- "Control"
s1$percent.mito <- PercentageFeatureSet(s1, pattern = "^mt-")
s1 <- RenameCells(s1, add.cell.id = 'Control')

# FAP_BiTE sample
s2_dir <- './data/SR006661-SR006239_TCell_FAP/count/sample_filtered_feature_bc_matrix/'
s2.data <- Read10X(data.dir = s2_dir)
s2 <- CreateSeuratObject(counts = s2.data)
s2$condition <- "FAP_BiTE"
s2$percent.mito <- PercentageFeatureSet(s2, pattern = "^mt-")
s2 <- RenameCells(s2, add.cell.id = 'FAPBiTE')
# s2 <- add_clonotype("./data/SR005215_CD45-2_Donor/vdj_t/", s2, "t")  # original commented

# Merge samples
sample <- merge(s1, y = c(s2))

# QC visualization & basic stats
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3, group.by = "condition")
sample

# Filters (features, mito, UMI)
sample <- subset(sample, subset = nFeature_RNA > 100 & nFeature_RNA < 10000 &
                                percent.mito < 10 & nCount_RNA < 50000)

# Post-filter QC
VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),
        ncol = 3, group.by = "condition")
sample

############################################################
## SCTransform, PCA/UMAP, neighbors, clustering
############################################################
DefaultAssay(sample) <- 'RNA'
sample <- SCTransform(sample, verbose = TRUE)
sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 100, verbose = TRUE)
sample <- RunUMAP(sample, reduction = "pca", dims = 1:30)
sample <- FindNeighbors(sample, reduction = "pca", dims = 1:30)
sample <- FindClusters(sample, graph.name = "SCT_snn", algorithm = 3,
                       resolution = c(0.1, 0.2, 0.3), verbose = FALSE)

# UMAP by cluster (res=0.3)
DimPlot(sample, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.3",
        cols = paletteDiscrete(unique(sample$SCT_snn_res.0.3), set = "stallion"),
        ncol = 2)

# DGE across res=0.3 clusters
Idents(sample) <- "SCT_snn_res.0.3"
DefaultAssay(sample) <- 'SCT'
rna.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./DE_SCT_snn_res.0.3.csv", quote = FALSE)

############################################################
## Focus on T Cells (clusters 3,4,6), recluster
############################################################
Idents(sample) <- "SCT_snn_res.0.3"
TCells <- subset(sample, idents = c("3","4","6"))

TCells <- RunUMAP(TCells, reduction = "pca", dims = 1:30)
TCells <- FindNeighbors(TCells, reduction = "pca", dims = 1:30)
TCells <- FindClusters(TCells, graph.name = "SCT_snn", algorithm = 3,
                       resolution = c(0.1, 0.2, 0.3, 0.4, 0.5), verbose = FALSE)

# UMAP by T-cell subclusters (res=0.4)
DimPlot(TCells, reduction = 'umap', label = FALSE, repel = TRUE, label.size = 2.5,
        group.by = "SCT_snn_res.0.4",
        cols = paletteDiscrete(unique(TCells$SCT_snn_res.0.4), set = "stallion"),
        ncol = 2)

# Composition of subclusters by condition
ggplot(TCells@meta.data,
       aes(x = condition, fill = SCT_snn_res.0.4)) +
  geom_bar(position = "fill") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(TCells$SCT_snn_res.0.4), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# DGE within T cell subclusters
Idents(TCells) <- "SCT_snn_res.0.4"
DefaultAssay(TCells) <- 'SCT'
rna.markers <- FindAllMarkers(TCells, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./Tcells/DE_SCT_snn_res.0.4.csv", quote = FALSE)

# DGE by condition in TCells
Idents(TCells) <- "condition"
DefaultAssay(TCells) <- 'SCT'
rna.markers <- FindAllMarkers(TCells, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.markers, file = "./Tcells/DE_condition.csv", quote = FALSE)

############################################################
## StarCat preprocessing (export counts + run external script)
############################################################
# (As in your original code; assumes convert_mouse_to_human_symbols() is available)
library(tidyverse)
library(data.table)
library(Matrix)
library(Seurat)
library(R.utils)
library(nichenetr)
library(SeuratDisk)

# Convert mouse symbols to human, rebuild object with same metadata
mouse_rna_matrix <- TCells@assays[["RNA"]]@counts
rownames(mouse_rna_matrix) <- mouse_rna_matrix %>% rownames() %>% convert_mouse_to_human_symbols()
mouse_rna_matrix <- mouse_rna_matrix %>% .[!is.na(rownames(mouse_rna_matrix)),
                                            !is.na(colnames(mouse_rna_matrix))]
mouse_coronary_new <- CreateSeuratObject(counts = mouse_rna_matrix)
mouse_coronary_new <- AddMetaData(object = mouse_coronary_new, metadata = TCells@meta.data)

# Write StarCat inputs (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
counts <- mouse_coronary_new@assays$RNA@counts
data_dir <- './Tcells/starCat/'

writeMM(counts, paste0(data_dir, 'matrix.mtx')); gzip(paste0(data_dir, 'matrix.mtx'))
barcodes <- colnames(counts)
readr::write_delim(as.data.frame(barcodes), paste0(data_dir, 'barcodes.tsv'),
                   col_names = FALSE); gzip(paste0(data_dir, 'barcodes.tsv'))
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names, type = "Gene Expression")
readr::write_delim(as.data.frame(features), paste0(data_dir, 'features.tsv'),
                   delim = "\t", col_names = FALSE); gzip(paste0(data_dir, 'features.tsv'))

# Build StarCat command
output_name <- 'tcells_bite_athero'
counts_fn <- paste0(data_dir, 'matrix.mtx.gz')
ref_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/In_vivo/Fap-targeting/Fap_BiTE_20wk/scRNAseq/TCR/Tcells/starCat/TCAT.V1/TCAT.V1.reference.tsv"
score_path <- "/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/Atherosclerosis/Projects/CITEseq_Atlas/In_vivo/Fap-targeting/Fap_BiTE_20wk/scRNAseq/TCR/Tcells/starCat/TCAT.V1/TCAT.V1.scores.yaml"

cmd <- paste0('/opt/anaconda3/bin/python ',
              '/opt/anaconda3/lib/python3.9/site-packages/starcat/starcat.py',
              ' --reference ', '"', ref_path, '"',
              ' --scores ', '"', score_path, '"',
              ' --counts ', '"', counts_fn, '"',
              ' --output-dir ', '"', data_dir, '"',
              ' --name ', '"', output_name, '"'
)
system(cmd)

# Load StarCat outputs
rf_usage_normalized <- read.table("./Tcells/starCat/tcells_bite_athero.rf_usage_normalized.txt")
scores <- read.table("./Tcells/starCat/tcells_bite_athero.scores.txt")

# Add score columns and RF usage as assay
TCells$starCat_ASA <- scores$ASA
TCells$starCat_Proliferation <- scores$Proliferation
TCells$starCat_ASA_binary <- scores$ASA_binary
TCells$starCat_Proliferation_binary <- scores$Proliferation_binary
TCells$starCat_Multinomial_Label <- scores$Multinomial_Label

rf_usage_normalized_assay <- CreateAssayObject(t(as.matrix(rf_usage_normalized)))
TCells@assays[["TCellRNorm"]] <- rf_usage_normalized_assay

# UMAPs colored by StarCat labels
DimPlot(TCells, label = FALSE, repel = FALSE, label.size = 3, reduction = 'umap',
        group.by = "starCat_Multinomial_Label",
        cols = paletteDiscrete(unique(TCells$starCat_Multinomial_Label), set = "stallion"))

# Average RF usage per T-cell subcluster
TCellRNorm_avg <- AverageExpression(TCells,
                                    features = rownames(TCells@assays[["TCellRNorm"]]),
                                    group.by = "SCT_snn_res.0.4")
pheatmap(t(TCellRNorm_avg$TCellRNorm), scale = "column", border_color = NA,
         col = paletteContinuous("solarExtra"), cluster_cols = TRUE)

# Binary flags to readable labels
TCells$starCat_ASA_binary[TCells$starCat_ASA_binary == "False"] <- "ASA_No"
TCells$starCat_ASA_binary[TCells$starCat_ASA_binary == "True"]  <- "ASA_Yes"
TCells$starCat_Proliferation_binary[TCells$starCat_Proliferation_binary == "False"] <- "Prolif_No"
TCells$starCat_Proliferation_binary[TCells$starCat_Proliferation_binary == "True"]  <- "Prolif_Yes"

# UMAPs of binary StarCat predictions
DimPlot(TCells, label = FALSE, repel = FALSE, label.size = 3, reduction = 'umap',
        group.by = "starCat_ASA_binary", cols = c("gray","red"))
DimPlot(TCells, label = FALSE, repel = FALSE, label.size = 3, reduction = 'umap',
        group.by = "starCat_Proliferation_binary", cols = c("gray","red"))

# Continuous StarCat scores as FeaturePlot
DefaultAssay(TCells) <- "RNA"
FeaturePlot(TCells, reduction = 'umap', features = "starCat_ASA") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,0.1))
FeaturePlot(TCells, reduction = 'umap', features = "starCat_Proliferation") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,0.1))

# Composition of ASA-positive per subcluster
ggplot(TCells@meta.data, aes(x = SCT_snn_res.0.4, fill = starCat_ASA_binary)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("gray","red")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Violin of ASA score by condition
VlnPlot(TCells, features = "starCat_ASA", group.by = "condition", ncol = 2)

############################################################
## T cell gene program scores (module scores)
############################################################
DefaultAssay(TCells) <- "SCT"

gene_modules <- list(
  Exhaustion = c("Pdcd1","Lag3","Ctla4","Tigit","Havcr2","Cd274","Tnfrsf9","Tnfrsf18",
                 "Tox","Nr4a1","Nr4a2","Nr4a3","Batf","Egr2","Irf4","Nfil3","Cebpb",
                 "Bhlhe40","Klf4","Runx1","Runx2","Stat3","Stat1","Ddit4","Hif1a",
                 "Tnfaip3","Cdkn1a","Ppp1r15a","Tsc22d3","Gadd45b","Ucp2","Hilpda"),
  Th17 = c("Il17a","Il17re","Rora","Stat3","Hif1a","Irf4","Cebpb"),
  Cytotoxic = c("Gzmk","Il2ra","Itgae","Ifngr2","Pik3cd","Ppp3ca","Cotl1","Il1r1",
                "Foxp1","Abi3","Tns4","Tespa1","Zbtb16"),
  Activation = c("Icos","Tnfrsf4","Cd40lg","Il21r","Grap2","Pik3ap1","Lax1","Tanc1",
                 "Tanc2","Nedd4","Il2rb","Rps28","Rps29"),
  Migration = c("Ccr2","Ccr8","Itgae","Itgav","Itpr1","Slc15a2","Cd6","Thy1","Itga1",
                "Rapgef4","Fyb","Sema4c"),
  Metabolic = c("Ddit4","Ucp2","Hilpda","Hif1a","Dusp1","Nfkbia","Tsc22d3","Pmaip1",
                "Stat5a","Epas1","Gadd45b","Tnfaip3"),
  Regulatory = c("Il1r2","Tnfaip3","Ctla4","Lag3","Pdcd1","Tigit","Tnfrsf18","Tnfsf11",
                 "Tsc22d3","Nfil3","Cebpb","Batf"),
  EarlyActivation = c("Junb","Fosl2","Dusp1","Dusp2","Dusp5","Ier5l","Btg1","Btg2","Nrgn"),
  Proliferation = c("Cdkn1a","Plk3","Cenpa","Runx1","Runx2","Stat3")
)

# Add module scores (suffix 1..9 appended automatically)
TCells <- AddModuleScore(TCells, features = gene_modules, name = names(gene_modules))

# Quick violins for each program by condition
VlnPlot(TCells, features = "Exhaustion1",      group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Th172",            group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Cytotoxic3",       group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Activation4",      group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Migration5",       group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Metabolic6",       group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Regulatory7",      group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "EarlyActivation8", group.by = "condition", ncol = 2)
VlnPlot(TCells, features = "Proliferation9",   group.by = "condition", ncol = 2)

# Save selected panels (pt.size=0 to declutter)
pdf("./Exhaustion.pdf", useDingbats = FALSE, width = 2, height = 3)
VlnPlot(TCells, features = "Exhaustion1", group.by = "condition", ncol = 2, pt.size = 0) + NoLegend()
dev.off()

pdf("./Th17.pdf", useDingbats = FALSE, width = 2, height = 3)
VlnPlot(TCells, features = "Th172", group.by = "condition", ncol = 2, pt.size = 0) + NoLegend()
dev.off()

pdf("./Cytotoxic.pdf", useDingbats = FALSE, width = 2, height = 3)
VlnPlot(TCells, features = "Cytotoxic3", group.by = "condition", ncol = 2, pt.size = 0) + NoLegend()
dev.off()

############################################################
## TCR: load contigs, combine, and repertoire visualizations
############################################################
# NOTE: The next few lines reference `tcells`, which is undefined earlier.
# If you meant `TCells`, uncomment the next line:
# tcells <- TCells

Idents(tcells) <- "condition"  # (will error if `tcells` not defined)
TCells_Control   <- subset(tcells, idents = "Control")
TCells_FAP_BiTE  <- subset(tcells, idents = "FAP_BiTE")

# Load filtered contigs and align barcodes to Seurat cell names
S1 <- read.csv("./data/SR006661-SR006239_TCell_Control/vdj_t/filtered_contig_annotations.csv")
S1$barcode2 <- paste("Control", S1$barcode, sep = "_")
S1 <- subset.data.frame(S1, S1$barcode2 %in% rownames(TCells_Control@meta.data))
S1$barcode2 <- NULL; S1

S2 <- read.csv("./data/SR006661-SR006239_TCell_FAP/vdj_t/filtered_contig_annotations.csv")
S2$barcode2 <- paste("FAPBiTE", S2$barcode, sep = "_")
S2 <- subset.data.frame(S2, S2$barcode2 %in% rownames(TCells_FAP_BiTE@meta.data))
S2$barcode2 <- NULL; S2

# Build scRepertoire list
contig_list <- list(S1, S2)
combined.TCR <- combineTCR(contig_list, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

# Repertoire summaries/plots
pdf("./clonalQuant.pdf", useDingbats = FALSE, width = 3, height = 3)
clonalQuant(combined.TCR, cloneCall = "strict", chain = "both", scale = TRUE)
dev.off()

clonalAbundance(combined.TCR, cloneCall = "gene", scale = FALSE)
clonalAbundance(combined.TCR, cloneCall = "gene", scale = TRUE)

pdf("./clonalLength.pdf", useDingbats = FALSE, width = 5, height = 3)
clonalLength(combined.TCR, cloneCall = "aa", chain = "both")
dev.off()

clonalCompare(combined.TCR, top.clones = 10, cloneCall = "aa", graph = "alluvial")

pdf("./clonalHomeostasis.pdf", useDingbats = FALSE, width = 5, height = 3)
clonalHomeostasis(combined.TCR, cloneCall = "gene")
dev.off()

pdf("./clonalHomeostasis.pdf", useDingbats = FALSE, width = 5, height = 3)
clonalHomeostasis(combined.TCR, cloneCall = "gene",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1,
                                Large = 0.3, Hyperexpanded = 1))
dev.off()

clonalProportion(combined.TCR, cloneCall = "gene")

pdf("./percentAA.pdf", useDingbats = FALSE, width = 5, height = 3)
percentAA(combined.TCR, chain = "TRB", aa.length = 20)
dev.off()

positionalEntropy(combined.TCR, chain = "TRB", aa.length = 20)

pdf("./percentGenes.pdf", useDingbats = FALSE, width = 5, height = 3)
percentGenes(combined.TCR, chain = "TRB", gene = "Vgene")
dev.off()

# PCA on V gene usage
df.genes <- percentGenes(combined.TCR, chain = "TRB", gene = "Vgene", exportTable = TRUE)
pc <- prcomp(df.genes)
df <- as.data.frame(cbind(pc$x[, 1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1); df$PC2 <- as.numeric(df$PC2)
ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = df[, 3]), shape = 21, size = 5) +
  guides(fill = guide_legend(title = "Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# PCA on VJ usage
percentVJ(combined.TCR, chain = "TRB")
df.genes <- percentVJ(combined.TCR, chain = "TRB", exportTable = TRUE)
pc <- prcomp(df.genes)
df <- as.data.frame(cbind(pc$x[, 1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1); df$PC2 <- as.numeric(df$PC2)
ggplot(df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = df[, 3]), shape = 21, size = 5) +
  guides(fill = guide_legend(title = "Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# K-mer content summaries
percentKmer(combined.TCR, cloneCall = "aa", chain = "TRB", motif.length = 3, top.motifs = 25)
percentKmer(combined.TCR, cloneCall = "nt", chain = "TRB", motif.length = 3, top.motifs = 25)

# Clonal diversity
pdf("./clonalDiversity.pdf", useDingbats = FALSE, width = 8, height = 3)
clonalDiversity(combined.TCR, cloneCall = "gene", n.boots = 20)
dev.off()

# Rarefaction using Shannon diversity (q = 1)
clonalRarefaction(combined.TCR, plot.type = 2, hill.numbers = 1, n.boots = 2)

TCells@meta.data

############################################################
## Map contig barcodes to Seurat cell names; combineExpression
############################################################
# Prefix barcodes to match `RenameCells` above
combined.TCR$S1$barcode <- paste0("Control_", combined.TCR$S1$barcode)
combined.TCR$S2$barcode <- paste0("FAPBiTE_", combined.TCR$S2$barcode)

# Subset TCells to those with TCR barcodes
TCells_2 <- subset(TCells, cells = c(combined.TCR$S1$barcode, combined.TCR$S2$barcode))

# Map condition to S1/S2 label
fun <- function(x) {
  if (x == "Control") {"S1"} else if (x == "FAP_BiTE") {"S2"}
}
TCells_2$condition2 <- mapply(fun, TCells_2$condition)

# Join TCR metrics back into Seurat
seurat <- combineExpression(combined.TCR, TCells_2, cloneCall = "gene", proportion = TRUE)
seurat@meta.data

# Heatmap of clone size vs subclusters
heatmap_matrix <- table(seurat$SCT_snn_res.0.4, seurat$cloneSize)
pheatmap(heatmap_matrix, color = viridis::viridis(100),
         main = "Clone Size Distribution Across T Cell States")

# Compositions
ggplot(seurat@meta.data, aes(x = cloneSize, fill = SCT_snn_res.0.4)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(seurat$SCT_snn_res.0.4), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggplot(seurat@meta.data, aes(x = condition, fill = cloneSize)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = as.vector(paletteDiscrete(unique(seurat$cloneSize), set = "stallion"))) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Exhaustion score vs clone size
pdf("./Exhaustion_violin.pdf", useDingbats = FALSE, width = 6, height = 4)
VlnPlot(seurat, features = "Exhaustion1", group.by = "cloneSize", ncol = 2, pt.size = 0)
dev.off()

# Save result
saveRDS(seurat, "./Tcells/Tcells_TCR.rds")

############################################################
## Additional panels (use `tcells` object; as in your code)
############################################################
# NOTE: the following use `tcells` which is not defined above; if needed:
# tcells <- TCells

pdf("./Tcell_dotplot.pdf", useDingbats = FALSE, width = 11, height = 3)
DotPlot(tcells,
        features = c("Cd8a","Gzmb","Ifng","Fasl","Tnfrsf9","Cd27","Cd28","Tigit","Hopx","Klrk1",
                     "Trdv4","Trdc","Il7r","Ltb4r1","Serpinb1a","Gpx1","Psap",
                     "Zbtb16","Il23r","Runx1",
                     "Tcrg-V5","Tcrg-V6","Junb","Cebpb","Crem","Icos","Klf4",
                     "Cd4","Cd2","Tcf7","Slamf6","Tnfrsf4","Tnfrsf18","Tnfsf8","Areg","Il1rl1",
                     "Bhlhe40","Rgs1","Rgs2","Rgs16","Ctla4","Thy1","Cd40lg","Nrgn",
                     "Sell","Ccr7","Bcl2","Lef1","Foxo1","Cd44"),
        group.by = "SCT_snn_res.0.4") + RotatedAxis()
dev.off()

pdf("./nFeature_RNA.pdf", useDingbats = FALSE, width = 4, height = 3)
VlnPlot(object = tcells, features = c("nFeature_RNA"),
        ncol = 4, group.by = "orig.ident", pt.size = 0)
dev.off()
