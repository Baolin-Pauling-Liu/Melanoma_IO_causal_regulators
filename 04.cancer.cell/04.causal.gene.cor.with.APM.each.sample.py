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
from scanpy.external.pp import dca as sc_dca
import os, gzip, shutil
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack, issparse

import scipy.sparse as sp # Add .A as a property if missing (CSR/CSC/COO just in case)
for cls in (sp.csr_matrix, sp.csc_matrix, sp.coo_matrix):
    if not hasattr(cls, "A"): cls.A = property(cls.toarray)

adata = sc.read_h5ad('/home/liubaoli/projects/01.melanoma/02.data/01.adata/10.cadata.h5ad')
adata.X = adata.raw.X
cancer_cell_index = pd.read_csv('/home/pauling/projects/01.melanoma/09.data.new/03.gene.module/01.malignant.cellid.csv')
adata_sub = adata[adata.obs.index.isin(cancer_cell_index['barcodekey'])]

adata_sub.layers["counts"] = adata_sub.X.copy()

vc = adata_sub.obs["Channel"].value_counts()
channels = vc[vc > 99].index.tolist()

#Protein-coding genes
gda = pd.read_csv('../02.data/02.gene/01.gene.csv')
set1 = set(adata_sub.var.index)
set2 = set(gda['gene_name'])
intersected_set = set1.intersection(set2)
intersected_list = list(intersected_set)
intersected_list = intersected_list

adata_sub = adata_sub[:, intersected_list].copy()
causal_gene = pd.read_csv('/home/pauling/projects/01.melanoma/09.data.new/19.causal.gene.new.Jan24.csv')

mhcI_genes  = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "B2M","TAP1","PSMB9","IRF1"]
mhcII_genes = ["HLA-DRA", "HLA-DRB1", "CD74", "HLA-DRB5",
               "HLA-DPA1", "HLA-DQA1", "HLA-DPB1", "HLA-DQB1", "HLA-DMA"]

APM_gene = mhcI_genes + mhcII_genes

BASE_DIR = "/home/liubaoli/projects/01.melanoma/02.data/11.gene.program/02.cancer.cell/57.APM"
JOBS = 59

for ch in channels:
    print(f"\n==> {ch}")
    outdir = os.path.join(BASE_DIR, ch)
    os.makedirs(outdir, exist_ok=True)

    # --- subset & QC ---
    adata_M = adata_sub[adata_sub.obs['Channel'] == ch].copy()
    sc.pp.filter_genes(adata_M, min_counts=3)
    sc.pp.filter_cells(adata_M, min_counts=3)
    
    # --- intersect gene universe with gda ---
    set1 = set(adata_M.var.index)
    set2 = set(gda['gene_name'])
    intersected_list = list(set1.intersection(set2))
    if len(intersected_list) < 5:
        print(f"[{ch}] too few genes after intersection with gda ({len(intersected_list)}). Skipping.")
        continue

    adata_M = adata_M[:, intersected_list]
    
    adata_M.var["mt"] = adata_M.var_names.str.startswith("MT-")
    # ribosomal genes
    adata_M.var["ribo"] = adata_M.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata_M.var["hb"] = adata_M.var_names.str.contains("^HB[^(P)]")
    
    sc.pp.calculate_qc_metrics(
        adata_M, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    
    # Normalizing to median total counts
    sc.pp.normalize_total(adata_M)
    # Logarithmize the data
    sc.pp.log1p(adata_M)
    
    # --- PCA / layers / extra filter ---
    sc.pp.pca(adata_M, n_comps=30)
    adata_M.layers["counts_csc"] = adata_M.layers["counts"].tocsc()
    if "n_counts" not in adata_M.obs:
        adata_M.obs["n_counts"] = np.asarray(adata_M.layers["counts_csc"].sum(axis=1)).ravel()
    sc.pp.filter_genes(adata_M, min_cells=5)

    # --- Hotspot core ---
    hs = hotspot.Hotspot(
        adata_M,
        layer_key="counts_csc",
        model='danb',
        latent_obsm_key="X_pca",
        umi_counts_obs_key="total_counts"
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30, approx_neighbors=False)
    hs.compute_autocorrelations(jobs=JOBS)

    # --- sets present in this sample ---
    apm_in = list(set(adata_M.var.index).intersection(set(APM_gene)))
    neg_in = list(set(adata_M.var.index).intersection(set(causal_gene['gene'])))

    # --- APM block (save NOW for this sample) ---
    if len(apm_in) >= 3:
        hs.compute_local_correlations(apm_in, jobs=JOBS)
        hs.create_modules(min_gene_threshold=4, core_only=False, fdr_threshold=0.05)
        apm_mod_path = os.path.join(outdir, f"{ch}.module.APM.csv")
        apm_cor_path = os.path.join(outdir, f"{ch}.APM.cor.csv")
        hs.modules.to_csv(apm_mod_path)
        hs.local_correlation_z.to_csv(apm_cor_path)
        print(f"[{ch}] wrote: {apm_mod_path}, {apm_cor_path}")
    else:
        print(f"[{ch}] APM genes in sample < 4 ({len(apm_in)}). Skipping APM module.")

    # --- negative set block (save NOW for this sample) ---
    if len(neg_in) >= 3:
        hs.compute_local_correlations(neg_in, jobs=JOBS)
        hs.create_modules(min_gene_threshold=2, core_only=False, fdr_threshold=0.05)
        neg_mod_path = os.path.join(outdir, f"{ch}.module.causal.gene.csv")
        neg_cor_path = os.path.join(outdir, f"{ch}.causal.cor.csv")
        hs.modules.to_csv(neg_mod_path)
        hs.local_correlation_z.to_csv(neg_cor_path)
        print(f"[{ch}] wrote: {neg_mod_path}, {neg_cor_path}")
    else:
        print(f"[{ch}] negative gene list in sample < 2 ({len(neg_in)}). Skipping negative module.")
