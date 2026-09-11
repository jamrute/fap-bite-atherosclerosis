#!/home/hk910097/miniconda3/envs/tangram-env/bin/python

import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import glob
import scanpy as sc
import scanpy.external as sce
import datetime
import re
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import tangram as tg
import matplotlib as mpl
import skimage
from anndata import AnnData
import pathlib
import tangram as tg
import math

## ------------------------- ##
## for giant
## vsm_dir = "/home/noco0013/ProcessedData/visium/new_out_with_reorient-image/CK*/outs/"
## for mi
vsm_dir = "/home/noco0013/ProcessedData/visium/Visium_*/outs/"
vsm_file = "filtered_feature_bc_matrix.h5"
rpt_file = "metrics_summary.csv"
samples = glob.glob(vsm_dir)
reports = glob.glob(vsm_dir+rpt_file)
samples.sort()
reports.sort()
## for giant
##outdir = "/home/noco0013/projects/hyojin-spatial/giant_cell/"
## for mi
outdir = "/home/noco0013/projects/hyojin-spatial/mi_visium/"


## ------------------------- ##
## scanpy 
## 	|_ https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html
## scanpy 
## 	|_ https://scanpy.readthedocs.io/en/latest/generated/scanpy.read_visium.html
## data integration 
## 	|_ https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html#Data-integration
## tangram
## 	|_ https://github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_without_squidpy.ipynb
## tangram x squidpy
## 	|_ https://squidpy.readthedocs.io/en/stable/external_tutorials/tutorial_tangram.html
## tangram code github
## 	|_ https://github.com/broadinstitute/Tangram/blob/master/tangram/plot_utils.py
## ------------------------- ##

## Load spatial data to AnnData by scanpy 
## ------------------------- ##
key = "Visium_"
adata_set =[]
qc1_set =[]
qc2_set =[]
sample_set =[]
for i in range(len(samples)):
	adata=sc.read_visium(samples[i], count_file=vsm_file)
	adata.var_names_make_unique()
	sid = [key+samples[i].split(key)[1].split("/")[0]]*adata.shape[0]
	adata.obs["sample_id"] = [key+samples[i].split(key)[1].split("/")[0]]*adata.shape[0]
	adata.var["mt"] = adata.var_names.str.startswith("MT-")
	sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
	adata_set.append(adata)
	qc_report = pd.read_csv(reports[i], sep=",")
	qc1_set.append( qc_report[["Number of Spots Under Tissue", "Valid UMIs", "Total Genes Detected" ]].values.tolist()[0] )
	qc2_set.append( adata.obs[["sample_id", "pct_counts_mt", "total_counts", "n_genes_by_counts", "log1p_total_counts","log1p_n_genes_by_counts"]].values.tolist())
	sample_set.append(sid[0])	


qc1_pd = pd.DataFrame(qc1_set, index=sample_set, columns = ["Number of Spots Under Tissue", "Valid UMIs", "Total Genes Detected" ])
qc2_pd = pd.DataFrame( [ sub2 for sub1 in qc2_set for sub2 in sub1 ], columns = ["sample", "pct_counts_mt", "total_counts", "n_genes_by_counts", "log1p_total_counts","log1p_n_genes_by_counts"])
qc1_pd["index"] = qc1_pd.index.tolist()
qc2_pd["index"] = qc2_pd.index.tolist()


## make bar plots
## ------------------------- ##
col1 = ["Number of Spots Under Tissue", "Valid UMIs", "Total Genes Detected" ]
fig, axs = plt.subplots(1, 3, figsize=(20, 4))
for n in range(len(col1)) :
	ax1=sns.barplot(data=qc1_pd, x='index', y=col1[n], dodge=True, ax=axs[n]);
	ax1.set_xticklabels(ax1.get_xticklabels(),rotation = 90)
	plt.savefig(os.path.join(outdir, "QC1.bar.across.samples.pdf"), bbox_inches='tight')
	#plt.tight_layout()

plt.close()


## make box plot
## ------------------------- ##
col2 = ["log1p_total_counts","log1p_n_genes_by_counts"]
fig, axs = plt.subplots(1, 2, figsize=(8, 4))
for n in range(len(col2)) :
	ax1 =sns.boxplot(x="sample", y=col2[n], data=qc2_pd, saturation=0.5, ax=axs[n])
	locs, labels = plt.xticks()
	ax1.set_xticklabels(ax1.get_xticklabels(),rotation = 90)
	plt.savefig(os.path.join(outdir, "QC2.box.log10.across.samples.pdf"), bbox_inches='tight')

plt.close()


col3 = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
fig, axs = plt.subplots(1, 3, figsize=(15, 4))
for n in range(len(col3)) :
        ax1= sns.boxplot(x="sample", y=col3[n], data=qc2_pd, saturation=0.5, ax=axs[n])
        ax1.set_xticklabels(ax1.get_xticklabels(),rotation = 90)
	plt.savefig(os.path.join(outdir, "QC2.box.across.samples.pdf"), bbox_inches='tight')

plt.close()


## make scatter plot
## ------------------------- ##
N = qc2_pd["sample"].value_counts().shape[0]
H = math.ceil(N/5)

wh =[]
for h in range(H):
	for w in range(5):
		wh.append([h,w])

## subplots(H,W)
## ------------------------- ##
fig, axs = plt.subplots(H, 5, figsize=(30, 5*H))
for n in range(len(sample_set)):
	dat = qc2_pd[qc2_pd["sample"]==sample_set[n]]
	sns.scatterplot(data=dat, x="total_counts", y="pct_counts_mt", ax=axs[wh[n][0], wh[n][1]]).set(title=sample_set[n])
	plt.savefig(os.path.join(outdir, "QC2.scatter.mt.by.count.across.samples.pdf"), bbox_inches='tight')

plt.close()

for n in range(len(sample_set)):
	dat = adata_set[n][adata_set[n].obs.sample_id == sample_set[n],:]  
	sc.pl.spatial(dat, library_id=sample_set[n], color = ["total_counts", "n_genes_by_counts",'pct_counts_mt'])
	plt.savefig(os.path.join(outdir, "QC3.spatial."+sample_set[n]+"dimplot.pdf"))
	plt.close()




## Filter and make qc plot
## ------------------------- ##
qc2_set_n =[]
cutoff_counts = [ 0, 1250, 0, 1000, 0, 
		1250, 500, 2000, 1000, 0,
		1250, 0, 1000, 1000, 1000, 
		500,  0,  0, 0, 0  
		 ]
for n in range(len(sample_set)):
	sc.pp.filter_cells(adata_set[n], min_counts=400)
	sc.pp.filter_cells(adata_set[n], min_counts=cutoff_counts[n])	
	sc.pp.filter_genes(adata_set[n], min_cells=10)
	qc2_set_n.append( adata_set[n].obs[["sample_id", "pct_counts_mt", "total_counts", "n_genes_by_counts", "log1p_total_counts","log1p_n_genes_by_counts"]].values.tolist())
	
qc2_pd_n = pd.DataFrame( [ sub2 for sub1 in qc2_set_n for sub2 in sub1 ], columns = ["sample", "pct_counts_mt", "total_counts", "n_genes_by_counts", "log1p_total_counts","log1p_n_genes_by_counts"])
qc2_pd_n["index"] = qc2_pd_n.index.tolist()

fig, axs = plt.subplots(H, 5, figsize=(30, 5*H))
for n in range(len(sample_set)):
        dat = qc2_pd_n[qc2_pd_n["sample"]==sample_set[n]]
        sns.scatterplot(data=dat, x="total_counts", y="pct_counts_mt", ax=axs[wh[n][0], wh[n][1]]).set(title=sample_set[n])
        plt.savefig(os.path.join(outdir, "post.QC2.scatter.mt.by.count.across.samples.pdf"), bbox_inches='tight')

plt.close()




## preprocess
## ------------------------- ##
for n in range(len(sample_set)):
	## Normalize, get HVG
	sc.pp.normalize_total(adata_set[n], inplace=True)
	sc.pp.log1p(adata_set[n])
	sc.pp.highly_variable_genes(adata_set[n], flavor="seurat", n_top_genes=2000)
	## PCA, UMAP 
	sc.pp.pca(adata_set[n])
	sc.pp.neighbors(adata_set[n])
	sc.tl.umap(adata_set[n])
	sc.tl.leiden(adata_set[n], resolution=1.5, key_added="clusters")
	## Visualize
	plt.rcParams["figure.figsize"] = (4, 4)
	sc.pl.umap(adata_set[n], color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
	plt.savefig(os.path.join(outdir, sample_set[n]+"_qc_after_filterout.pdf"))
	plt.close()
	## Visualization in spatial coordinates
	plt.rcParams["figure.figsize"] = (8, 8)
	sc.pl.spatial(adata_set[n], img_key="hires", color=["total_counts", "n_genes_by_counts"])
	plt.savefig(os.path.join(outdir, sample_set[n]+".sp.pdf"))
	plt.close()
	sc.pl.spatial(adata_set[n], img_key="hires", color="clusters", size=1.5)
	plt.savefig(os.path.join(outdir, sample_set[n]+".sp.umap.pdf"))
	plt.close()






## ref data 
## ------------------------------------------ ##
## R1_f = "/home/noco0013/projects/mi_hearts/snRNA/snRNA-seq-submission.h5ad"
## R1 = scvi.data.read_h5ad(R1_f)
## R1.obs['batch'] = R1.obs['patient'].astype(str) + '_' + R1.obs['seq_batch'].astype(str)
## R1.obs.batch.value_counts()
## R1.obs.cell_type_original.value_counts()
## R1.write_h5ad(filename=R1_f+".cellmap.h5ad")
## ------------------------------------------ ##
new_f = ## PLEASE ADD YOUR REF DATA ##"/home/noco0013/projects/cellmap/data/heart/internal_mi/snRNA-seq-submission.cellmap.h5ad"
ad_sc = sc.read(new_f)
ad_sc = ad_sc
ad_sc.raw = ad_sc
celltype_column_name = 'cell_type_original'


## for a quick test or block
## ad_sc = ad_sc[:,0:100]
## ------------------------------------------ ##
## preprocess
sc.pp.normalize_total(ad_sc, inplace=True)
sc.pp.log1p(ad_sc)
sc.pp.highly_variable_genes(ad_sc, flavor="seurat", n_top_genes=2000)
sc.tl.rank_genes_groups(ad_sc, groupby=celltype_column_name, use_raw=False)


## use the intersection of the highly variable genes.
markers_df = pd.DataFrame(ad_sc.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
genes_sc = np.unique(markers_df.melt().value.values)


## Tangram
## Deconvolution and mapping
for i in range(len(samples)):
	## ------------------------- ##
	## tangram set up 
	## ------------------------- ##
	ad_sp = adata_set[i]
	ad_sp.obs['sample'] = list(ad_sp.uns['spatial'].keys())[0]
	##
	## find mitochondria-encoded (MT) genes
	## ------------------------- ##
	ad_sp.var['MT_gene'] = [gene.startswith('MT-') for gene in ad_sp.var.index]
	##
	## remove MT genes for spatial mapping (keeping their counts in the object)
	## ------------------------- ##
	ad_sp.obsm['MT'] = ad_sp[:, ad_sp.var['MT_gene'].values].X.toarray()
	ad_sp = ad_sp[:, ~ad_sp.var['MT_gene'].values]
	## ------------------------- ##
	## run Tangram
	## 	|_ 1.  to run Tangram at cell level
	## 	|_ 2.  to run Tangram at cluster level
	##
	## 1. run Tangram at cell level
	## ------------------------- ##
	## pp_adatas finds the common genes between adata_sc, adata_sp, 
	## and saves them in two adatas.uns for mapping and analysis later.
	## ------------------------- ##
	## If genes=None, Tangram maps using all genes shared by the two datasets
	## ------------------------- ##
	genes_st = ad_sp.var_names.values
	genes = list(set(genes_sc).intersection(set(genes_st)))
	tg.pp_adatas(ad_sc, ad_sp, genes=None)
	## ------------------------- ##
	## ad_map, is a cell-by-voxel structure 
	## where ad_map.X[i, j] gives the probability for cell i to be in voxel j
	## ------------------------- ##
	## plot_cell_annotation_sc parameter (1), (2), (3)
	## (1) Spot Size & (2) Scale Factor : shoule be None when ad_sp.uns['spatial'] exist
	## (3) ax : A matplotlib axes object. Only works if plotting a single component. 
	## ------------------------- ##
	## def plot_cell_annotation_sc(
	##    adata_sp, 
	##    annotation_list, 
	##    x="x", 
	##    y="y", 
	##    spot_size=None, 
	##    scale_factor=None, 
	##    perc=0,
	##    ax=None
	##):
	## ------------------------- ##
	## output saved in 
	## df = ad_sp.obsm["tangram_ct_pred"][annotation_list]
	## ------------------------- ##
	ad_map = tg.map_cells_to_space(ad_sc, ad_sp)
	## ------------------------- ##
	tg.project_cell_annotations(ad_map, ad_sp, annotation=celltype_column_name)
	annotation_list = list(pd.unique(ad_sc.obs[celltype_column_name]))
	tg.plot_cell_annotation_sc(ad_sp, annotation_list,x='x', y='y',spot_size= None, perc=0.1)
	plt.savefig(os.path.join(outdir, ad_sp.obs['sample'][0]+".tangram.pdf"), bbox_inches='tight')
	plt.close()
	## ad_ge is a voxel-by-gene AnnData, similar to spatial data ad_sp, 
	## but where gene expression has been projected from the single cells
	## ------------------------- ##
	ad_ge = tg.project_genes(ad_map, ad_sc)
	## 
	##
	tg.plot_training_scores(ad_map, bins=10, alpha=.5)
	plt.savefig(os.path.join(outdir, ad_sp.obs["sample"][0]+".tangram.plot_training_scores.pdf"), bbox_inches='tight')
	plt.close()
	## ------------------------- ##
	## write
	## ------------------------- ##
	ad_sp.write_h5ad(filename=outdir+ad_sp.obs["sample"][0]+".tangram.h5ad")
	##
	## 
	##tg.plot_genes_sc(genes, adata_measured=ad_sp, adata_predicted=ad_ge, spot_size=50, perc = 0.001, return_figure=False)
	##plt.savefig(os.path.join(outdir, ad_sp.obs["sample"][0]+".tangram.plot_genes_sc.pdf"), bbox_inches='tight')
	##plt.close()
	##
	##
	##df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_sp, ad_sc)
	##df_all_genes.to_csv( ad_sp.obs["sample"][0]+".tangram.compare_spatial_geneexp.txt", sep="\t")
	##tg.plot_auc(df_all_genes)
	##plt.savefig(os.path.join(outdir, ad_sp.obs['sample'][0]+".tangram.plot_auc.pdf"), bbox_inches='tight')
	##plt.close()
	## 
		








