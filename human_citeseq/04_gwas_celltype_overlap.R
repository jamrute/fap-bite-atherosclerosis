############################################################
# GWAS–DGE overlap analysis (deduped libraries + comments)
# - Keeps your original logic
# - Removes duplicate/unused library calls
# - Adds clear, step-by-step comments
############################################################

## ---- Libraries (unique + sufficient) ----
# tidyverse gives you dplyr, ggplot2, readr, tibble, etc.
library(tidyverse)

## ==========================================================
## 1) Load DGE table and keep significant genes
## ==========================================================
# Expecting columns like: gene, p_val_adj, cluster, etc.
DGE <- readr::read_delim("./DE_SCT_cell.type.csv", delim = ",",
                         escape_double = FALSE, trim_ws = TRUE)
DGE

# Filter to FDR-significant hits and pull unique gene list
DGE <- dplyr::filter(DGE, p_val_adj < 0.05)
diff_exp_genes <- unique(DGE$gene)
length(diff_exp_genes)

## ==========================================================
## 2) Load GWAS gene list
## ==========================================================
# One column named "gene" in cad_gwas_genes.txt
GWAS <- readr::read_delim("./cad_gwas_genes.txt", delim = ",",
                          escape_double = FALSE, trim_ws = TRUE)
gwas_genes <- GWAS$gene
gwas_genes

## ==========================================================
## 3) Intersect and focus DGE on GWAS-overlapping genes
## ==========================================================
intersecting_genes <- base::Reduce(intersect, list(diff_exp_genes, gwas_genes))
length(intersecting_genes)

DGE_subset <- dplyr::filter(DGE, gene %in% intersecting_genes)
DGE_subset

## ==========================================================
## 4) Count overlaps by cluster and visualize
## ==========================================================
# How many GWAS-overlap genes per cluster?
cad_gwas_overlap <- data.frame(sort(table(DGE_subset$cluster), decreasing = TRUE))
cad_gwas_overlap

# Lollipop-style plot (as in your code)
pdf("./cad_gwas_overlap.pdf", useDingbats = FALSE, width = 3.5, height = 2.5)

# Fix cluster order for the x-axis to your specified biology-driven order
cad_gwas_overlap$Var1 <- factor(
  cad_gwas_overlap$Var1,
  levels = c("BCells","TCells","pDC","Mast","Proliferating","Fibroblast1",
             "Myeloid","Glia","Lymphatic","Fibroblast2","SMCPericyte",
             "Endothelium","ModSMC")
)

ggplot(cad_gwas_overlap, aes(x = Var1, y = Freq)) +
  geom_segment(aes(xend = Var1, yend = 0.5)) +
  geom_point(size = 2, color = "black") +
  coord_flip() +
  theme_bw() +
  xlab("")

dev.off()
