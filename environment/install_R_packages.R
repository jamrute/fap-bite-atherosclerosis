################################################################################
# Install the R packages used in this repository.
# Paper : Amrute et al., Science (2026) | doi:10.1126/science.adx1736
#
# Usage : Rscript environment/install_R_packages.R
# Notes : Scripts use Seurat v4-style assay access (e.g. obj@assays$RNA@counts,
#         slot = "data"); with Seurat v5, set
#         options(Seurat.object.assay.version = "v3") before creating objects.
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran <- c(
  # single-cell core
  "Seurat", "Signac", "sctransform", "harmony",
  # data handling
  "dplyr", "tidyr", "tibble", "readr", "tidyverse", "data.table", "Matrix",
  "reshape2", "R.utils",
  # plotting
  "ggplot2", "ggpubr", "ggrepel", "ggsci", "patchwork", "cowplot", "pheatmap",
  "viridis", "RColorBrewer", "scales", "corrplot", "UpSetR", "gt", "clustree",
  "scCustomize",
  # analysis utilities
  "igraph", "uwot", "cluster", "compositions", "msigdbr", "RcppML"
)

bioc <- c(
  "clusterProfiler", "DOSE", "enrichplot", "ReactomePA",
  "org.Hs.eg.db", "org.Mm.eg.db", "biomaRt",
  "progeny", "dorothea", "viper", "UCell", "fgsea", "scran",
  "Nebulosa", "scRepertoire"
)

github <- c(
  "GreenleafLab/ArchR",            # paletteDiscrete() / paletteContinuous()
  "mojaveazure/seurat-disk",       # SeuratDisk: h5Seurat <-> h5ad
  "carmonalab/GeneNMF",            # NMF meta-programs
  "rpolicastro/scProportionTest",  # cell-type proportion permutation test
  "saeyslab/nichenetr",            # convert_mouse_to_human_symbols()
  "neurorestore/Augur"             # cell-type perturbation prioritisation
)

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
for (repo in github) remotes::install_github(repo, upgrade = "never")

message("Done. Record versions with sessionInfo() for reproducibility.")
