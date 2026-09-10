############################################################
# SMC/Pericyte/Fibroblast analysis + Owens mapping + ATAC label transfer
# - Libraries deduplicated & completed (SeuratDisk, org.Hs.eg.db, progeny, etc.)
# - Original flow preserved with light commenting
# - Uses ArchR palettes: paletteDiscrete / paletteContinuous
############################################################

## ---- Libraries (unique & sufficient) ----
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(harmony)
library(ArchR)         # paletteDiscrete/paletteContinuous
library(ggsci)
library(SeuratDisk)
library(readr)
library(pheatmap)
library(viridis)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)
library(progeny)
library(Nebulosa)
library(dorothea)
library(tibble)
library(tidyr)
library(viper)
library(ggrepel)
library(scProportionTest)

############################################################
# Load SMC/Pericyte/Fibroblast object and basic prep
############################################################
sample <- readRDS("SMCPericyte_Fibroblast.rds")

DefaultAssay(sample) <- "RNA"
sample <- NormalizeData(sample)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
sample <- ScaleData(sample)

DefaultAssay(sample) <- "ADT"
sample <- NormalizeData(sample, normalization.method = 'CLR', margin = 2) %>% ScaleData()

# Clustering (RNA graph assumed)
sample <- FindClusters(sample, graph.name = "RNA_snn", algorithm = 3,
                       resolution = c(0.4,0.5,0.6,0.7), verbose = FALSE)

saveRDS(sample, "./SMCPericyte_Fibroblast.rds")

# UMAP labels (expects rna.umap present in object)
DimPlot(sample, reduction = 'rna.umap', group.by = 'RNA_snn_res.0.7',
        label.size = 4, label = TRUE,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.7), set = "stallion"))

# Annotate states
fun <- function(x) {
  if (x == "0") "Fib2"
  else if (x == "1") "SMC2"
  else if (x == "2") "SMC1"
  else if (x == "3") "Fib5"
  else if (x == "4") "SMC3"
  else if (x == "5") "CMC"
  else if (x == "6") "FMC"
  else if (x == "7") "Fib1"
  else if (x == "8") "Pericyte"
  else if (x == "9") "Fib6"
  else if (x == "10") "Fib4"
  else if (x == "11") "Fib7"
  else if (x == "12") "SMC4"
  else if (x == "13") "Fib3"
}
sample$cell.state <- mapply(fun, sample$RNA_snn_res.0.7)
sample$cell.state <- factor(sample$cell.state,
  levels = c("Pericyte","SMC1","SMC2","SMC3","SMC4","FMC","CMC","Fib1","Fib2","Fib3","Fib4","Fib5","Fib6","Fib7"))

DimPlot(sample, reduction = 'rna.umap', group.by = 'cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(sample$cell.state), set = "stallion"))

# DE (RNA)
DefaultAssay(sample) <- 'RNA'
Idents(sample) <- "cell.state"
rna.rnamarkers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(rna.rnamarkers, file ="./DE_RNA_cell.state.csv", quote = FALSE)

# Top markers heatmap (RNA)
rna.rnamarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sample, features = top10$gene, assay = "RNA",
          group.colors = paletteDiscrete(unique(sample$cell.state), set = "stallion"),
          size = 2, angle = 90) + NoLegend() +
  scale_fill_gradientn(colours = paletteContinuous("solarExtra"))
ggsave(filename = "heatmap_DE_RNA_cell.state.pdf")

# DE (ADT)
DefaultAssay(sample) <- 'ADT'
Idents(sample) <- "cell.state"
adt.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(adt.markers, file ="./DE_ADT_cell.state.csv", quote = FALSE)

adt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sample, features = top10$gene, assay = "ADT",
          group.colors = paletteDiscrete(unique(sample$cell.state), set = "stallion"),
          size = 2, angle = 90) + NoLegend() +
  scale_fill_gradientn(colours = paletteContinuous("blueYellow"))
ggsave(filename = "heatmap_ADT.pdf")

saveRDS(sample, "SMCPericyte_Fibroblast.rds")

# Quick markers/plots in annotated object
sample <- readRDS("smc_fib_annotated.rds")
DefaultAssay(sample) <- "SCT"
DotPlot(sample, features = c("IFIH1","FAP"), group.by = "cell.state")
FeaturePlot(sample, features = "ITGA2")

############################################################
# Owens mouse mapping to SMC reference
############################################################
owens_global_mapped <- readRDS("/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/CAD_Project/Projects/CITE-seq Atlas/analysis/final_analysis/Owens_Mapping/reference_mapped.rds")
mouse_old <- readRDS("/Users/jamrute/Library/CloudStorage/Box-Box/Macbook_Files/Grad_School/Primary_Projects/CAD_Project/Projects/CITE-seq Atlas/analysis/final_analysis/Owens_Mapping/mouse_humanHomolog_normalized.rds")
mouse_coronary <- UpdateSeuratObject(mouse_old)

mouse_coronary_new <- CreateSeuratObject(counts = mouse_coronary@assays[["RNA"]]@counts)
mouse_coronary_new <- AddMetaData(mouse_coronary_new, metadata = mouse_coronary@meta.data)
mouse_coronary_new$predicted.celltype <- owens_global_mapped$predicted.celltype

Idents(mouse_coronary_new) <- "predicted.celltype"
owens_SMC <- subset(mouse_coronary_new, idents = c("ModSMC","SMCPericyte","Fibroblast1"))

DefaultAssay(owens_SMC) <- "RNA"
owens_SMC <- NormalizeData(owens_SMC)
owens_SMC <- FindVariableFeatures(owens_SMC, selection.method = "vst", nfeatures = 3000)
owens_SMC <- ScaleData(owens_SMC)

# Reference mapping from annotated human SMC object
sample <- RunSPCA(sample, assay = 'RNA', graph = 'RNA_snn')
anchors <- FindTransferAnchors(reference = sample, query = owens_SMC,
                               normalization.method = "LogNormalize",
                               reference.reduction = "spca",
                               dims = 1:50, reference.assay = "RNA")
owens_SMC <- MapQuery(anchorset = anchors, query = owens_SMC, reference = sample,
                      refdata = list(celltype = "RNA_snn_res.0.3", predicted_ADT = "ADT"),
                      reference.reduction = "spca", reduction.model = "rna.umap")

DimPlot(owens_SMC, reduction = "ref.umap", group.by = "predicted.celltype", label.size = 4,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.3), set = "stallion"), label = FALSE)

Idents(owens_SMC) <- "orig.cell"
owens_SMC$orig.ident <- factor(owens_SMC$orig.ident, levels = c("SMC","Endothelial","Other"))
DimPlot(owens_SMC, reduction = "ref.umap", group.by = "orig.cell",
        label = FALSE, label.size = 3, repel = TRUE, pt.size = 0.01, shuffle = TRUE) +
  scale_colour_manual(values = c("grey","grey","red"), na.value = "grey")
DimPlot(owens_SMC, reduction = "ref.umap", group.by = "orig.cell",
        label = FALSE, label.size = 3, repel = TRUE, pt.size = 0.01, shuffle = TRUE) +
  scale_colour_manual(values = c("blue","grey","grey"), na.value = "grey")

owens_SMC$SEURATPROJECT <- "SMCs"
ggplot(owens_SMC@meta.data, aes(x = SEURATPROJECT, fill = predicted.celltype)) +
  geom_bar(position = "fill") + theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("#D51F26","#272E6A","#208A42","#89288F","#8A9FD1","#C06CAB")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

table(owens_SMC$predicted.celltype, owens_SMC$orig.cell)

# h5ad export
owens_SMC$orig.cell <- as.character(owens_SMC$orig.cell)
SaveH5Seurat(owens_SMC, filename = "./owens_mapping/SMC_ref_mapped.h5Seurat")
Convert("./owens_mapping/SMC_ref_mapped.h5Seurat", dest = "h5ad")
saveRDS(owens_SMC, "./owens_mapping/SMC_ref_mapped.rds")

# Prediction score heatmap
DefaultAssay(owens_SMC) <- "prediction.score.celltype"
owens_SMC@assays[["prediction.score.celltype"]]@counts <- owens_SMC@assays[["prediction.score.celltype"]]@data
Idents(owens_SMC) <- "predicted.celltype"
avg <- AverageExpression(owens_SMC,
                         features = rownames(owens_SMC@assays[["prediction.score.celltype"]]),
                         assays = "prediction.score.celltype",
                         group.by = "predicted.celltype")
avg[[1]] <- avg[[1]][c("0","1","2","3","4","5","6","7"),]
mat <- as.matrix(avg)[[1]]
pheatmap(mat, scale = "none", col = viridis(240), cexCol = 0.5, cellwidth = 20,
         cluster_rows = FALSE, fontsize_row = 6, fontsize_col = 6,
         cluster_cols = FALSE, legend = TRUE, border_color = NA, cellheight = 20,
         breaks = seq(0,1,by=1/240))

############################################################
# Pathway analysis across clusters (GO)
############################################################
d <- read.csv("./v2_DE_RNA_snn_res.0.3.csv") %>%
     filter(avg_log2FC > 0.58, p_val_adj < 0.05) %>%
     select(gene, cluster)

eg <- bitr(as.character(d$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
d2 <- d %>% filter(gene %in% eg$SYMBOL)
d_ent <- merge(d2, eg, by.x = "gene", by.y = "SYMBOL") %>% select(ENTREZID, cluster)
geneList <- unstack(d_ent)
ck <- compareCluster(geneCluster = geneList, fun = enrichGO, OrgDb = "org.Hs.eg.db")
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
dotplot(ck, font.size = 8) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

############################################################
# PROGENy pathway activities (by cell.state)
############################################################
Idents(sample) <- "cell.state"
CellsClusters <- data.frame(Cell = names(Idents(sample)),
                            CellType = as.character(Idents(sample)),
                            stringsAsFactors = FALSE)
sample <- progeny(sample, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE)
sample <- Seurat::ScaleData(sample, assay = "progeny")

progeny_scores_df <- as.data.frame(t(GetAssayData(sample, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) %>%
  inner_join(CellsClusters)

summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  select(-std) %>% spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue","white","red"))(paletteLength)
progenyBreaks <- c(seq(min(summarized_progeny_scores_df), 0, length.out = ceiling(paletteLength/2) + 1),
                   seq(max(summarized_progeny_scores_df)/paletteLength,
                       max(summarized_progeny_scores_df), length.out = floor(paletteLength/2)))
pheatmap(t(summarized_progeny_scores_df[,-1]), fontsize = 14, fontsize_row = 10,
         color = myColor, breaks = progenyBreaks, main = "PROGENy (500)",
         angle_col = 45, border_color = NA, scale = "row")

############################################################
# Palantir I/O + FDL embedding
############################################################
write.csv(as.matrix(sample[["RNA"]]@scale.data),
          file = "./palantir/smc_fibro_RNA_normalized.txt", quote = FALSE)
write.csv(sample@meta.data, file = "./palantir/smc_fibro_meta.csv", quote = TRUE)

sample_meta <- read.csv2('./palantir/Stroma_palantir_meta_data.csv', header = TRUE, sep = ',', row.names = 1)
sample_2 <- AddMetaData(sample, sample_meta)
sample_2@meta.data$pseudotime <- as.numeric(as.character(sample_2@meta.data$pseudotime))
sample_2@meta.data$entropy   <- as.numeric(as.character(sample_2@meta.data$entropy))

fdl <- read.csv2('./palantir/Stroma_fdl.csv', header = TRUE, sep = ',', row.names = 1)
fdl$x <- as.double(fdl$x); fdl$y <- as.double(fdl$y)
colnames(fdl) <- paste0("FDL_", 1:2)
sample_2[["fdl"]] <- CreateDimReducObject(embeddings = as.matrix(fdl), key = "FDL_")

DimPlot(sample_2, reduction = 'fdl', group.by = 'RNA_snn_res.0.7', label.size = 4,
        label = FALSE, cols = paletteDiscrete(unique(sample_2$RNA_snn_res.0.7), set = "stallion"))

FeaturePlot(sample_2, reduction = "fdl", features = c("pseudotime")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,0.75))
FeaturePlot(sample_2, reduction = "fdl", features = c("entropy")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,0.5))

sample_meta$pseudotime <- as.numeric(sample_meta$pseudotime)
sample_meta$entropy    <- as.numeric(sample_meta$entropy)
sample_meta$X4 <- as.numeric(sample_meta$X4)
sample_meta$X5 <- as.numeric(sample_meta$X5)
sample_meta$X1 <- as.numeric(sample_meta$X1)
sample_meta$ClusterName <- as.character(sample$RNA_snn_res.0.7)
sample_2 <- AddMetaData(sample, sample_meta)

FeaturePlot(sample_2, reduction = "fdl", features = c("X4")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,1))
FeaturePlot(sample_2, reduction = "fdl", features = c("X5")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,1))
FeaturePlot(sample_2, reduction = "fdl", features = c("X1")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,1))

ggplot(sample_meta, aes(x = pseudotime, y = entropy, color = ClusterName)) +
  geom_point(size = 1, alpha = 0.1) +
  scale_color_manual(values = paletteDiscrete(unique(sample$RNA_snn_res.0.7), set = "stallion")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

FeaturePlot(sample_2, reduction = "rna.umap", features = c("pseudotime")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,0.75))
FeaturePlot(sample_2, reduction = "rna.umap", features = c("entropy")) +
  scale_color_gradientn(colors = paletteContinuous("horizon"), oob = scales::squish, limits = c(0,0.5))

VlnPlot(sample, features = "MFGE8",
        cols = paletteDiscrete(unique(sample_2$RNA_snn_res.0.3), set = "stallion"),
        sort = TRUE, pt.size = 0) + NoLegend()
VlnPlot(sample_2, features = "entropy",
        cols = paletteDiscrete(unique(sample_2$RNA_snn_res.0.3), set = "stallion"),
        sort = TRUE, pt.size = 0) + NoLegend()

############################################################
# Label Transfer ATAC (Miller) to RNA states
############################################################
atac_mapped_SMC <- readRDS("./atac_mapping/miller_atac_SMC_labelTransfer.rds")
DimPlot(atac_mapped_SMC, reduction = "umap.atac", group.by = "predicted.id", label.size = 4,
        cols = paletteDiscrete(unique(sample$RNA_snn_res.0.3), set = "stallion"), label = FALSE)

# Summaries of prediction scores
df <- atac_mapped_SMC@meta.data[, c("prediction.score.5","prediction.score.0","prediction.score.1",
                                    "prediction.score.4","prediction.score.2","prediction.score.3",
                                    "prediction.score.6","prediction.score.7","predicted.id")]
df2 <- data.frame(sort(unique(df$predicted.id)))
df2$"5" <- tapply(df$prediction.score.5, list(df$predicted.id), mean)
df2$"0" <- tapply(df$prediction.score.0, list(df$predicted.id), mean)
df2$"1" <- tapply(df$prediction.score.1, list(df$predicted.id), mean)
df2$"4" <- tapply(df$prediction.score.4, list(df$predicted.id), mean)
df2$"2" <- tapply(df$prediction.score.2, list(df$predicted.id), mean)
df2$"3" <- tapply(df$prediction.score.3, list(df$predicted.id), mean)
df2$"6" <- tapply(df$prediction.score.6, list(df$predicted.id), mean)
df2$"7" <- tapply(df$prediction.score.7, list(df$predicted.id), mean)

df3 <- df2[,-1]; rownames(df3) <- df2[,1]
df3 <- df3[, c("0","1","2","3","4","5","6","7")]
unique(atac_mapped_SMC$predicted.id)
df4 <- t(df3)[, c("0","1","2","3","4","5","6")]

pheatmap(df4, color = viridis(250), border_color = NA, cellwidth = 20, scale = "none",
         cluster_rows = FALSE, cluster_cols = FALSE, cellheight = 20)

df <- atac_mapped_SMC@meta.data[, c("prediction.score.max","predicted.id")]
df2 <- data.frame(sort(unique(df$predicted.id)))
df2$maxScore <- tapply(df$prediction.score.max, list(df$predicted.id), mean)
df3 <- df2[,-1]; rownames(df3) <- df2[,1]
df4 <- t(df3)[, c("0","1","2","3","4","5","6")]; df4["7"] <- 0
df5 <- t(df4)[, c("0","1","2","3","4","5","6","7")]

pheatmap(df5, color = viridis(250), border_color = NA, cellwidth = 20, scale = "none",
         cellheight = 20, cluster_rows = FALSE, cluster_cols = FALSE)

FeaturePlot(atac_mapped_SMC, features = c("prediction.score.0"), reduction = "umap.atac",
            cols = c("lightgrey","darkred")) & theme(plot.title = element_text(size = 10))
FeaturePlot(atac_mapped_SMC, features = c("prediction.score.3"), reduction = "umap.atac",
            cols = c("lightgrey","darkred")) & theme(plot.title = element_text(size = 10))
FeaturePlot(atac_mapped_SMC, features = c("prediction.score.1"), reduction = "umap.atac",
            cols = c("lightgrey","darkred")) & theme(plot.title = element_text(size = 10))

# Peak-level DE across predicted.id
DefaultAssay(atac_mapped_SMC) <- 'Peaks'
Idents(atac_mapped_SMC) <- "predicted.id"
de.peaks <- FindAllMarkers(atac_mapped_SMC, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.25)
open_c3 <- rownames(de.peaks[de.peaks$avg_log2FC > 0.5, ])
closest_genes_c3 <- ClosestFeature(atac_mapped_SMC, regions = open_c3)

############################################################
# TF enrichment (DoRothEA + VIPER)
############################################################
Idents(sample) <- "RNA_snn_res.0.3"
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>% filter(confidence %in% c("A","B","C"))
sample <- run_viper(sample, regulon,
                    options = list(method = "scale", minsize = 4, eset.filter = FALSE,
                                   cores = 1, verbose = FALSE))
DefaultAssay(sample) <- "dorothea"
sample <- ScaleData(sample)

Idents(sample) <- "RNA_snn_res.0.3"
viper_scores_df <- GetAssayData(sample, slot = "scale.data", assay = "dorothea") %>%
  data.frame(check.names = FALSE) %>% t()
CellsClusters <- data.frame(cell = names(Idents(sample)),
                            cell_type = as.character(Idents(sample)),
                            check.names = FALSE)
viper_scores_clusters <- viper_scores_df %>%
  data.frame() %>% rownames_to_column("cell") %>% gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>% summarise(avg = mean(activity), std = sd(activity))
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>% mutate(var = var(avg)) %>%
  ungroup() %>% top_n(1000, var) %>% distinct(tf)
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  select(-std) %>% spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE)

############################################################
# Gene module scores & signatures
############################################################
DefaultAssay(sample) <- "RNA"
# SMC signature
expdata <- GetAssayData(sample)
Pop1 <- c("ACTA2","CNN1","MYL9","TPM2","MYH11","TAGLN","SOST","PPP1R14A","COL18A1","ITIH4")
zz <- which(tolower(rownames(expdata)) %in% tolower(Pop1))
geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp))); geneExp[is.nan(geneExp)] <- 0
sample@meta.data$SMC <- colSums(geneExp) / length(zz)
FeaturePlot(sample, features = "SMC", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,1))
plot_density(sample, features = "RUNX1", reduction = 'rna.umap')

# Chondrocyte-like score
Pop1 <- c("OMD","LUM","HAPLN1","COMP","FMOD","FAP","ITGA10","COL8A1","RUNX1","LTBP2","ENPP1")
zz <- which(tolower(rownames(expdata)) %in% tolower(Pop1))
geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp))); geneExp[is.nan(geneExp)] <- 0
sample@meta.data$chondro <- colSums(geneExp) / length(zz)
FeaturePlot(sample, features = "chondro", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,1))

DefaultAssay(sample) <- "SCT"
FeaturePlot(sample, reduction = "rna.umap", features = "FAP") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"), oob = scales::squish, limits = c(0,1))

# Example ADT feature on a separate object 'cor'
cor <- readRDS("smc_fib_annotated.rds")
DefaultAssay(cor) <- "ADT"
FeaturePlot(cor, reduction = "rna.umap", features = "FAP.1") +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,2))

# Inflammation/foam-like signature (IPCOX)
DefaultAssay(sample) <- "RNA"
expdata <- GetAssayData(sample)
Pop1 <- c("MMP12","IL1B","PTGS2","EREG","CCL5","IL1A","CXCL3","CXCL8","EDNRB","STC1",
          "HMOX1","AKR1C2","MMP1","DUSP6","NAMPT","CXCL2","AKR1C1","SERPINB2","FADS1",
          "THBD","ADAMTS4","PLIN2","MT2A","LIF","AKR1B1")
zz <- which(tolower(rownames(expdata)) %in% tolower(Pop1))
geneExp <- as.matrix(expdata[zz, ])
geneExp <- t(scale(t(geneExp))); geneExp[is.nan(geneExp)] <- 0
sample@meta.data$IPCOX <- colSums(geneExp) / length(zz)
FeaturePlot(sample, features = "IPCOX", reduction = 'rna.umap') +
  scale_color_gradientn(colors = c("blue","turquoise2","yellow","red","red4"),
                        oob = scales::squish, limits = c(0,0.5))

############################################################
# Proportion tests by sex on annotated object
############################################################
cor <- readRDS("smc_fib_annotated.rds")
prop_test <- sc_utils(cor)
prop_test <- permutation_test(prop_test, cluster_identity = "cell.state",
                              sample_1 = "F", sample_2 = "M",
                              sample_identity = "Sex")
permutation_plot(prop_test)

############################################################
# RUNX & other plots
############################################################
Idents(sample) <- "cell.state"
SMC_lineage <- subset(sample, idents = c("SMC1","SMC2","SMC3","SMC4","FMC","CMC"))
DotPlot(sample, features = c("TNFRSF11B","FAP"), group.by = "cell.state") + RotatedAxis()
DefaultAssay(cor) <- "RNA"
VlnPlot(cor, features = "IGFBP7", group.by = "cell.state", pt.size = 0)

############################################################
# Bulk DEG signatures to module scores
############################################################
g1 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_FoamCell1d_vs_BM1d.csv", sep = ",")
g2 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_FoamCell8d_vs_BM8d.csv", sep = ",")
g3 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_FoamCell8d_vs_FoamCell1d.csv", sep = ",")
g4 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_oxLDL1d_vs_BM1d.csv", sep = ",")
g5 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_oxLDL8d_vs_BM8d.csv", sep = ",")
g6 <- read.csv2("./Fw_ Bulk RNA-seq project_ SMC in foam cell media/DEG_oxLDL8d_vs_oxLDL1d.csv", sep = ",")

DefaultAssay(sample) <- "RNA"
sample <- AddModuleScore(sample, list(g1$Symbol), name = 'FoamCell1d_vs_BM1d')
sample <- AddModuleScore(sample, list(g2$Symbol), name = 'FoamCell8d_vs_BM8d')
sample <- AddModuleScore(sample, list(g3$Symbol), name = 'FoamCell8d_vs_FoamCell1d')
sample <- AddModuleScore(sample, list(g4$Symbol), name = 'oxLDL1d_vs_BM1d')
sample <- AddModuleScore(sample, list(g5$Symbol), name = 'oxLDL8d_vs_BM8d')
sample <- AddModuleScore(sample, list(g6$Symbol), name = 'oxLDL8d_vs_oxLDL1d')

Idents(sample) <- "cell.state"
sample2 <- subset(sample, idents = c("Pericyte","SMC1","SMC2","SMC3","SMC4","CMC","FMC"))
DotPlot(sample2,
        features = c("FoamCell1d_vs_BM1d1","FoamCell8d_vs_BM8d1","FoamCell8d_vs_FoamCell1d1",
                     "oxLDL1d_vs_BM1d1","oxLDL8d_vs_BM8d1","oxLDL8d_vs_oxLDL1d1"),
        group.by = "cell.state") + RotatedAxis()

sample2 <- AddModuleScore(sample2, list(intersect(g2$Symbol, g5$Symbol)),
                          name = 'oxLDL8d_AND_FoamCell8d')
DotPlot(sample2, features = "RUNX1", group.by = "cell.state") + RotatedAxis()

DefaultAssay(sample) <- "SCT"
DotPlot(sample, features = c("RUNX2","MFGE8","TP53","RB1","CDKN1A","CDKN2A")) + RotatedAxis()

############################################################
# VSMC compartment w/ Palantir metadata
############################################################
sample_meta <- read.csv2('./palantir/Stroma_palantir_meta_data.csv', header = TRUE, sep = ',', row.names = 1)
sample2 <- AddMetaData(sample, sample_meta)
sample2@meta.data$pseudotime <- as.numeric(as.character(sample2@meta.data$pseudotime))
sample2@meta.data$entropy   <- as.numeric(as.character(sample2@meta.data$entropy))

Idents(sample2) <- "cell.state"
sample2 <- subset(sample2, idents = c("Pericyte","SMC1","SMC2","SMC3","SMC4","FMC","CMC"))
sample2 <- RunUMAP(sample2, reduction = "harmony_rna", dims = 1:50,
                   reduction.name = 'rna.umap2', reduction.key = 'rnaUMAP_', assay = "RNA")

DimPlot(sample2, reduction = 'rna.umap', group.by = 'cell.state',
        label.size = 4, label = FALSE,
        cols = paletteDiscrete(unique(sample2$cell.state), set = "stallion"))

FeaturePlot(sample2, reduction = "rna.umap", features = "pseudotime") +
  scale_color_gradientn(colors = paletteContinuous("horizon"),
                        oob = scales::squish, limits = c(0,0.5))
FeaturePlot(sample2, reduction = "rna.umap", features = "entropy") +
  scale_color_gradientn(colors = paletteContinuous("horizon"),
                        oob = scales::squish, limits = c(0,0.5))

DefaultAssay(sample2) <- "SCT"
FeaturePlot(sample2, reduction = "rna.umap", features = "FAP") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,1))
FeaturePlot(sample2, reduction = "rna.umap", features = "MYH11") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,3))

############################################################
# FAP protein–RNA correlation analysis (SMC lineage only)
############################################################
find_fap_correlations <- function(seurat_obj, genes = NULL, min_correlation = 0.3, p_value_cutoff = 0.05) {
  expr_matrix <- GetAssayData(seurat_obj[["RNA"]], slot = "data")
  adt_matrix  <- GetAssayData(seurat_obj[["ADT"]], slot = "data")
  fap_expr    <- adt_matrix["FAP.1", ]

  if (!is.null(genes)) expr_matrix <- expr_matrix[rownames(expr_matrix) %in% genes, , drop = FALSE]
  expr_matrix <- expr_matrix[rownames(expr_matrix) != "FAP", , drop = FALSE]

  cors <- cor(t(as.matrix(expr_matrix)), as.matrix(fap_expr), method = "spearman")
  n <- ncol(expr_matrix)
  t_stat <- cors * sqrt((n - 2) / (1 - cors^2))
  p_values <- 2 * pt(-abs(t_stat), df = n - 2)

  results <- data.frame(gene = rownames(expr_matrix),
                        correlation = as.numeric(cors),
                        p_value = as.numeric(p_values),
                        row.names = NULL)
  results$p_adj <- p.adjust(results$p_value, method = "BH")

  sig <- results[abs(results$correlation) >= min_correlation & results$p_adj <= p_value_cutoff, ]
  sig[order(abs(sig$correlation), decreasing = TRUE), ]
}

plot_fap_correlations <- function(correlation_results, top_n = 30) {
  plot_data <- head(correlation_results, top_n)
  ggplot(plot_data, aes(x = reorder(gene, correlation), y = correlation)) +
    geom_bar(stat = "identity",
             fill = ifelse(plot_data$correlation > 0, "#4C78A8", "#E45756")) +
    theme_minimal() + coord_flip() +
    labs(x = "Gene", y = "Correlation with FAP Protein",
         title = paste("Top", top_n, "genes correlating with FAP expression")) +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 11))
}

Idents(sample) <- "cell.state"
VSMC <- subset(sample, idents = c("Pericyte","SMC1","SMC2","SMC3","SMC4","FMC","CMC"))
fap_correlations <- find_fap_correlations(VSMC)
head(fap_correlations); fap_correlations

pos_cors <- subset(fap_correlations, correlation > 0)
neg_cors <- subset(fap_correlations, correlation < 0)
plot_fap_correlations(fap_correlations)
write.csv(fap_correlations, file = "./FAP_protein_correlations_protein.csv", quote = FALSE)

DefaultAssay(sample) <- "ADT"
FeaturePlot(sample, reduction = "rna.umap", features = "F3.1") +
  scale_color_gradientn(colors = paletteContinuous("blueYellow"),
                        oob = scales::squish, limits = c(0,4))

############################################################
# Marker dotplots & THY1 RNA/ADT
############################################################
sample <- readRDS("smc_fib_annotated.rds")
d <- read.csv("./DE_RNA_cell.state.csv") %>% filter(avg_log2FC > 0.58, p_val_adj < 0.05)
d %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top10; top10

DefaultAssay(sample) <- "SCT"
pdf("./marker_dotplot.pdf", useDingbats = FALSE, width = 16, height = 4)
DotPlot(sample, features = unique(top10$gene), group.by = "cell.state",
        col.min = 0, cols = c("lightgrey","red")) + RotatedAxis()
dev.off()

fap_correlations <- read.csv2("./FAP_protein_correlations_protein.csv", sep = ",")
fap_correlations$correlation <- as.double(fap_correlations$correlation)
pdf("./FAP_protein_protein_corr.pdf", useDingbats = FALSE, width = 2, height = 2)
plot_fap_correlations(fap_correlations)
dev.off()

DefaultAssay(sample) <- "SCT"
FeaturePlot(sample, reduction = "rna.umap", features = "LUM") +
  scale_color_gradientn(colors = paletteContinuous("solarExtra"),
                        oob = scales::squish, limits = c(0,3))

sample <- readRDS("smc_fib_annotated.rds")
Idents(sample) <- "cell.state"
sample <- subset(sample, idents = c("Pericyte","SMC1","SMC2","SMC3","SMC4","FMC","CMC"))

DefaultAssay(sample) <- "SCT"
DotPlot(sample, features = c("LOXL1","PALLD","MYO9B","NAV1","SBF2","FURIN","AMOTL2")) + RotatedAxis()

Idents(sample) <- "cell.state"
sample <- subset(sample, idents = c("Pericyte","SMC1","SMC2","SMC3","SMC4","FMC","CMC"))

pdf("./THY1_Protein.pdf", useDingbats = FALSE, width = 3.5, height = 4)
DefaultAssay(sample) <- "ADT"
DotPlot(sample, features = "THY1.1") + RotatedAxis()
dev.off()

DefaultAssay(sample) <- "SCT"
DotPlot(sample, features = c("CD68","SPP1","LGALS3","PLIN2","MSR1","CD36","OLR1","TREM2","APOE",
                             "LAPTM5","CD74","GRN","FABP5","AIF1","LAMP2","ABCA1")) + RotatedAxis()
FeaturePlot(sample, features = "NRPC3")

sample <- readRDS("smc_fib_annotated.rds")
pdf("./THY1_RNA.pdf", useDingbats = FALSE, width = 3.5, height = 4)
DefaultAssay(sample) <- "SCT"
DotPlot(sample, features = "THY1", group.by = "cell.state") + RotatedAxis()
dev.off()

pdf("./THY1_protein.pdf", useDingbats = FALSE, width = 3.5, height = 4)
DefaultAssay(sample) <- "ADT"
DotPlot(sample, features = "THY1.1", group.by = "cell.state") + RotatedAxis()
dev.off()

############################################################
# Collapse to QSMC vs modSMC; RNA/ADT DE and intersections
############################################################
fun <- function(x) {
  if (x == "Pericyte") "Pericyte"
  else if (x %in% c("SMC1","SMC2","SMC3","SMC4")) "QSMC"
  else if (x %in% c("FMC","CMC")) "modSMC"
  else if (grepl("^Fib", x)) "Fib"
}
sample$cell.state2 <- mapply(fun, sample$cell.state)
Idents(sample) <- "cell.state2"
SMCs <- subset(sample, idents = c("QSMC","modSMC"))

DefaultAssay(SMCs) <- 'RNA'
Idents(SMCs) <- "cell.state2"
rna.rnamarkers <- FindAllMarkers(SMCs, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
write.csv(rna.rnamarkers, file ="./DE_RNA_Protein/DE_RNA_cell.state.csv", quote = FALSE)

DefaultAssay(SMCs) <- 'ADT'
Idents(SMCs) <- "cell.state2"
prot.rnamarkers <- FindAllMarkers(SMCs, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
write.csv(prot.rnamarkers, file ="./DE_RNA_Protein/DE_ADT_cell.state.csv", quote = FALSE)

# Volcano-like prep (flip sign so positive = modSMC)
de <- prot.rnamarkers
de$avg_log2FC <- ifelse(de$cluster == "QSMC", -de$avg_log2FC, de$avg_log2FC)
de$diffexpressed <- "NO"
de$diffexpressed[de$avg_log2FC >  0.58 & de$p_val_adj < 0.05] <- "UP"
de$diffexpressed[de$avg_log2FC < -0.58 & de$p_val_adj < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]

pdf("./volcano_protein.pdf", useDingbats = FALSE, width = 5, height = 4)
ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed)) +
  geom_point() + theme_minimal() +
  scale_color_manual(values = c("black","grey","black"))
dev.off()

# Intersections QSMC/modSMC across RNA vs ADT
rna.q <- subset(rna.rnamarkers, p_val_adj < 0.05 & avg_log2FC > 0.25 & cluster == "QSMC")
adt.q <- subset(prot.rnamarkers, p_val_adj < 0.05 & avg_log2FC > 0.25 & cluster == "QSMC")
intersect(rna.q$gene, sub("\\.1$","", adt.q$gene))

rna.m <- subset(rna.rnamarkers, p_val_adj < 0.05 & avg_log2FC > 0.25 & cluster == "modSMC")
adt.m <- subset(prot.rnamarkers, p_val_adj < 0.05 & avg_log2FC > 0.25 & cluster == "modSMC")
intersect(rna.m$gene, sub("\\.1$","", adt.m$gene))

# Merge per-modSMC genes (RNA vs ADT)
rna <- read.csv("./DE_RNA_Protein/DE_RNA_cell.state.csv") %>% subset(p_val_adj < 0.05)
rna.mod <- subset(rna, cluster == "modSMC")[, c("gene","avg_log2FC")]

protein <- read.csv("./DE_RNA_Protein/DE_ADT_cell.state.csv")
protein$gene <- sub("\\.1$","", protein$gene)
protein <- subset(protein, p_val_adj < 0.05)
prot.mod <- subset(protein, cluster == "modSMC")[, c("gene","avg_log2FC")]

result <- rna.mod %>% inner_join(prot.mod, by = "gene", suffix = c(".df1",".df2"))
write.csv(result, file = "./DE_RNA_Protein/merged.csv", quote = FALSE)
