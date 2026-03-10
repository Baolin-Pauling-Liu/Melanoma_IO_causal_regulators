import numpy as np
import pandas as pd
import scanpy as sc 
import anndata as an
import bbknn
import seaborn as sns
import scanpy.external as sce

adata = sc.read_h5ad("../02.data/01.adata/adata.pro.clustering.h5ad")
adata.X = adata.raw.X   #Raw counts data

gda = pd.read_csv('../02.data/02.gene/01.gene.csv')  # Protein-coding genes
set1 = set(adata.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)

adata = adata[:,intersected_list]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, 'Sample')  #Batch correction

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1)

adata.write_h5ad("../02.data/01.adata/adata.pro.clustering.h5ad")
