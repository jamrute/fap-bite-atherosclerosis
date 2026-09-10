############################################
# Libraries
library(Seurat)        # core single-cell workflow
library(SeuratDisk)    # Convert/Load h5ad <-> h5Seurat
library(ggplot2)       # plotting
library(patchwork)     # plot composition
library(dplyr)         # data wrangling
library(readr)         # read_tsv for metadata
library(remotes)       # optional: install Github pkgs (GeneNMF)
# remotes::install_github("carmonalab/GeneNMF")  # install if not available
library(GeneNMF)       # NMF-based program discovery
library(scCustomize)   # convenience wrappers; used for SCTransform call here
library(Matrix)        # dgCMatrix class support
library(RcppML)        # ML backends used by GeneNMF
library(reshape2)      # melt() for ggplot heatmap
library(pheatmap)      # heatmap for Jaccard matrices

## ==========================================================
## 1) Convert AnnData -> h5Seurat and load Seurat object
## ==========================================================
options(future.globals.maxSize = 50000 * 1024^2)  # allow up to ~50 GB for futures

# Convert h5ad to h5Seurat (overwrites if exists)
Convert("new_SIKANDER.h5ad", dest = "h5seurat", overwrite = TRUE)

# Load as Seurat object
merged <- LoadH5Seurat("new_SIKANDER.h5seurat")

## ==========================================================
## 2) Read and attach metadata
##    - keep only barcodes present in 'merged'
##    - ensure uniqueness by barcode
## ==========================================================
meta <- read_tsv("metadata.tsv", show_col_types = FALSE) |> as.data.frame()
meta <- meta |> dplyr::filter(barcode %in% colnames(merged)) |> dplyr::distinct(barcode, .keep_all = TRUE)
rownames(meta) <- meta$barcode
meta$barcode <- NULL
merged <- AddMetaData(merged, metadata = meta[colnames(merged), , drop = FALSE])

## ==========================================================
## 3) SCTransform normalization prior to program discovery
## ==========================================================
merged <- SCTransform(merged, assay = "RNA", verbose = TRUE)

## Alias (optional convenience, kept from your code)
seu <- merged

## ==========================================================
## 4) NMF decomposition and multi-sample metaprograms
## ==========================================================
ndim <- 30  # rank for runNMF

# Find variable features then run NMF on SCT assay
seu <- FindVariableFeatures(seu, nfeatures = 3000)
seu <- runNMF(seu, k = ndim, assay = "SCT")

# Split by sample and run multiNMF across a range of k
seu.list <- SplitObject(seu, split.by = "sample")
geneNMF.programs <- multiNMF(seu.list, assay = "SCT", slot = "data", k = 4:9, nfeatures = 3000)

# Derive metaprograms (consensus across runs)
geneNMF.metaprograms <- getMetaPrograms(
  geneNMF.programs,
  nMP = 10,               # number of metaprograms
  weight.explained = 0.7, # variance explained threshold
  max.genes = 100         # max genes per MP
)

# Quick visualization & inspection
plotMetaPrograms(geneNMF.metaprograms)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)

## ==========================================================
## 5) Collect MP gene lists and save
## ==========================================================
MP1_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP1);  MP1_genes$MP  <- "MP1";  colnames(MP1_genes)  <- c("gene","MP")
MP2_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP2);  MP2_genes$MP  <- "MP2";  colnames(MP2_genes)  <- c("gene","MP")
MP3_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP3);  MP3_genes$MP  <- "MP3";  colnames(MP3_genes)  <- c("gene","MP")
MP4_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP4);  MP4_genes$MP  <- "MP4";  colnames(MP4_genes)  <- c("gene","MP")
MP5_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP5);  MP5_genes$MP  <- "MP5";  colnames(MP5_genes)  <- c("gene","MP")
MP6_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP6);  MP6_genes$MP  <- "MP6";  colnames(MP6_genes)  <- c("gene","MP")
MP7_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP7);  MP7_genes$MP  <- "MP7";  colnames(MP7_genes)  <- c("gene","MP")
MP8_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP8);  MP8_genes$MP  <- "MP8";  colnames(MP8_genes)  <- c("gene","MP")
MP9_genes  <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP9);  MP9_genes$MP  <- "MP9";  colnames(MP9_genes)  <- c("gene","MP")
MP10_genes <- data.frame(geneNMF.metaprograms$metaprograms.genes$MP10); MP10_genes$MP <- "MP10"; colnames(MP10_genes) <- c("gene","MP")

MP_genes <- dplyr::bind_rows(MP1_genes, MP2_genes, MP3_genes, MP4_genes,
                             MP5_genes, MP6_genes, MP7_genes, MP8_genes,
                             MP9_genes, MP10_genes)

write.csv(MP_genes, file = "./Sikander_MP_genes.csv", quote = FALSE)

## ==========================================================
## 6) Jaccard overlap analysis: CA vs Sikander MPs
##    - builds pairwise and cross-dataset heatmaps
## ==========================================================
# Read MP gene tables (update paths if needed)
ca_data       <- read.csv("MP_genes_CA.csv",       stringsAsFactors = FALSE)
sikander_data <- read.csv("Sikander_MP_genes.csv", stringsAsFactors = FALSE)

# Label source dataset and combine
ca_data$dataset        <- "CA"
sikander_data$dataset  <- "Sikander"
combined_data          <- rbind(ca_data, sikander_data)

# Unique metaprogram identifiers per dataset
combined_data$MP_unique <- paste(combined_data$dataset, combined_data$MP, sep = "_")

# Jaccard index helper
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

# Per-MP gene lists
unique_mps <- unique(combined_data$MP_unique)
mp_gene_lists <- lapply(unique_mps, function(mp) combined_data$gene[combined_data$MP_unique == mp])
names(mp_gene_lists) <- unique_mps

# Pairwise Jaccard matrix (all MPs across both datasets)
n_mps <- length(unique_mps)
jaccard_matrix <- matrix(0, nrow = n_mps, ncol = n_mps,
                         dimnames = list(unique_mps, unique_mps))
for (i in 1:n_mps) {
  for (j in 1:n_mps) {
    jaccard_matrix[i, j] <- jaccard_index(mp_gene_lists[[i]], mp_gene_lists[[j]])
  }
}

# Heatmap (pheatmap)
pheatmap(jaccard_matrix,
         main = "Jaccard Index Heatmap: Gene Overlap Between MPs",
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 8, cellwidth = 30, cellheight = 30, angle_col = 45)

# Alternative ggplot heatmap
jaccard_df <- reshape2::melt(jaccard_matrix)
colnames(jaccard_df) <- c("MP1", "MP2", "Jaccard_Index")
ggplot_heatmap <- ggplot(jaccard_df, aes(x = MP2, y = MP1, fill = Jaccard_Index)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                       midpoint = 0.5, name = "Jaccard\nIndex") +
  geom_text(aes(label = round(Jaccard_Index, 2)), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "Jaccard Index Heatmap: Gene Overlap Between MPs",
       x = "Molecular Phenotype",
       y = "Molecular Phenotype")
print(ggplot_heatmap)

# Summary stats & MP gene counts
cat("Summary of Jaccard Index values:\n")
print(summary(as.vector(jaccard_matrix[upper.tri(jaccard_matrix)])))
cat("\nNumber of genes per MP:\n")
gene_counts <- sapply(mp_gene_lists, length)
print(gene_counts)

# Cross-dataset Jaccard (CA MPs vs Sikander MPs only)
ca_mps        <- unique(ca_data$MP)
sikander_mps  <- unique(sikander_data$MP)
cross_jaccard_matrix <- matrix(0, nrow = length(ca_mps), ncol = length(sikander_mps),
                               dimnames = list(paste("CA", ca_mps, sep = "_"),
                                               paste("Sikander", sikander_mps, sep = "_")))
for (i in 1:length(ca_mps)) {
  ca_genes <- ca_data$gene[ca_data$MP == ca_mps[i]]
  for (j in 1:length(sikander_mps)) {
    sikander_genes <- sikander_data$gene[sikander_data$MP == sikander_mps[j]]
    cross_jaccard_matrix[i, j] <- jaccard_index(ca_genes, sikander_genes)
  }
}

# Save cross-dataset heatmap and report maxima
pdf("./Cross-Dataset Jaccard Index.pdf", useDingbats = FALSE, width = 5, height = 5)
pheatmap(cross_jaccard_matrix,
         main = "Cross-Dataset Jaccard Index: CA vs Sikander MPs",
         color = colorRampPalette(c("white", "lightcoral", "red", "darkred"))(100),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = FALSE, cellwidth = 10, cellheight = 10,
         angle_col = 45, border_color = NA)
dev.off()

cat("\nCross-dataset Jaccard Index summary (CA vs Sikander):\n")
print(summary(as.vector(cross_jaccard_matrix)))

# Highest cross-dataset overlap(s)
max_jaccard  <- max(cross_jaccard_matrix)
max_indices  <- which(cross_jaccard_matrix == max_jaccard, arr.ind = TRUE)
cat("\nHighest cross-dataset overlap:\n")
for (k in 1:nrow(max_indices)) {
  i <- max_indices[k, 1]; j <- max_indices[k, 2]
  cat(sprintf("%s vs %s: Jaccard Index = %.3f\n",
              rownames(cross_jaccard_matrix)[i],
              colnames(cross_jaccard_matrix)[j],
              max_jaccard))
}

# End of script
