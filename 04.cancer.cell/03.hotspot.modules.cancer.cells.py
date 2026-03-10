import pandas as pd
import scanpy as sc 
import anndata as an
import bbknn
import seaborn as sns
import scanpy.external as sce
import hotspot
import sys
import numpy as np
import pickle
from hotspot.plots import local_correlation_plot
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage

adata = sc.read_h5ad('../02.data/01.adata/10.cadata.harmony.h5ad')
adata.X = adata.raw.X
cancer_cell_index = pd.read_csv('/home/pauling/projects/01.melanoma/09.data.new/03.gene.module/01.malignant.cellid.csv')
adata_sub = adata[adata.obs.index.isin(cancer_cell_index['barcodekey'])]

#Protein-coding genes
gda = pd.read_csv('../02.data/02.gene/01.gene.csv')
set1 = set(adata.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)

adata_sub = adata_sub[:, intersected_list]

adata_sub.layers["counts"] = adata_sub.X.copy()

sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)

adata_sub.layers["counts_csc"] = adata_sub.layers["counts"].tocsc()
sc.pp.filter_genes(adata_sub, min_cells=50)
hs = hotspot.Hotspot(
    adata_sub,
    layer_key="counts_csc",
    model='danb',
    latent_obsm_key="X_pca_harmony",
    umi_counts_obs_key="n_counts"
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=20)

hs_results = hs.compute_autocorrelations(jobs=50)
sc.pp.highly_variable_genes(adata_sub, min_mean=0.05, max_mean=5, min_disp=0.1)

adata_sub.var['highly_variable'].value_counts()

hvg_genes = adata_sub.var_names[adata_sub.var['highly_variable']]
lcz = hs.compute_local_correlations(hvg_genes, jobs=50)

hs.local_correlation_z.to_csv("/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/50.gene.cor.csv")

modules = hs.create_modules(
    min_gene_threshold=10, core_only=False, fdr_threshold=0.05
)
module_scores = hs.calculate_module_scores()

hs.plot_local_correlations(vmin = -50, vmax = 50, yticklabels = False) #gene module plot -- heatmap

hs.modules.to_csv("/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/50.HVG.module.csv")
module_scores.to_csv("/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/50.HVG.module.score.csv")
