import pandas as pd 
import scanpy as sc 
import anndata as an
import bbknn
import numpy as np
import seaborn as sns
import scanpy.external as sce
import hotspot
import sys
import scipy.sparse as sp

adata = sc.read_h5ad("../02.data/01.adata/11.bcell.h5ad")  # Replace with corresponding cell types

gda = pd.read_csv('../02.data/02.gene/01.gene.csv')
set1 = set(adata.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)

gda = pd.read_csv('../02.data/02.gene/01.gene.csv')
set1 = set(adata.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)

adata_sub = adata[:, intersected_list]

assert sp.issparse(adata_sub.X), "adata.X must be sparse to use this method efficiently."

# Step 1: Create a sparse binary matrix (1 if expression > 0)
binary_expr = adata_sub.X.copy()
binary_expr.data = np.ones_like(binary_expr.data)

# Step 2: Get the sample labels
samples = adata_sub.obs['Sample'].values

# Step 3: Unique samples
unique_samples = np.unique(samples)

# Step 4: Prepare a result container
result = {}

# Step 5: For each sample, calculate fraction
for sample in unique_samples:
    idx = np.where(samples == sample)[0]  # indices of cells from this sample
    if len(idx) == 0:
        continue
    
    # Subset the binary matrix
    sub_binary_expr = binary_expr[idx, :]
    
    # Sum across cells (axis=0), and divide by number of cells
    fraction = sub_binary_expr.sum(axis=0).A1 / len(idx)  # .A1 flattens to 1D array
    
    result[sample] = fraction

# Step 6: Convert to DataFrame
fraction_expressing = pd.DataFrame(result, index=adata_sub.var_names).T

cell_counts = adata_sub.obs['Sample'].value_counts()

# Step 2: Select samples where count > 9
samples_with_enough_cells = cell_counts[cell_counts > 9].index

medians = fraction_expressing.loc[samples_with_enough_cells,:].quantile(0.85, axis=0)

# Step 3: Filter genes with median > 0.04
filtered_genes = medians[medians > 0.04].index

adata_sub = adata_sub[:, filtered_genes]
raw = adata_sub.raw[:,adata_sub.var.index]
adata_sub.X = raw.X

cda_filt = adata_sub

cda_filt.layers["counts"] = cda_filt.X.copy()
sc.pp.normalize_total(cda_filt, target_sum=1e4)
sc.pp.log1p(cda_filt)
cda_filt.layers["log_normalized"] = cda_filt.X.copy()

cda_filt.layers["counts_csc"] = cda_filt.layers["counts"].tocsc()
hs = hotspot.Hotspot(
    cda_filt,
    layer_key="counts_csc",
    model='danb',
    latent_obsm_key="X_pca",
    umi_counts_obs_key="n_counts"
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

lcz = hs.compute_local_correlations(filtered_genes, jobs=50)

modules = hs.create_modules(
    min_gene_threshold=10, core_only=False, fdr_threshold=0.05
)

hs.plot_local_correlations(vmin = -40, vmax = 40)

module_scores = hs.calculate_module_scores()

module_scores.to_csv("/home/liubaoli/projects/01.melanoma/02.data/08.hotspot/04.B.cell/Bcell.module.score.csv")
lcz.to_csv("/home/liubaoli/projects/01.melanoma/02.data/08.hotspot/04.B.cell/Bcell.gene.cor.csv")
hs.modules.to_csv("/home/liubaoli/projects/01.melanoma/02.data/08.hotspot/04.B.cell/Bcell.gene.module.csv")
