import numpy as np
import pandas as pd
import scanpy as sc 
import anndata as an
import bbknn
import seaborn as sns
import scanpy.external as sce
import scrublet as scr

adata = sc.read_h5ad("../02.data/01.adata/adata.pro.clustering.h5ad")
adata.X = adata.raw.X.copy()

samples = pd.read_csv("/home/pauling/projects/01.melanoma/09.data.new/02.sample.info/sample.clinial.meta.csv")
adata_sub = adata[adata.obs['Sample'].isin(samples['sample'])]

gda = pd.read_csv('../02.data/02.gene/01.gene.csv')
set1 = set(adata_sub.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)

adata_sub = adata_sub[:, intersected_list]
count_expr = adata_sub.X

scrub = scr.Scrublet(count_expr)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata_sub.obs['doublet'] = predicted_doublets
adata_sub.obs['doublet_score'] = doublet_scores

adata_sub.obs.to_csv("/home/liubaoli/projects/01.melanoma/02.data/01.adata/01.round1.new.clustering.csv")
