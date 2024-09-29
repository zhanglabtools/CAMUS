# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 20:52:37 2024

@author: Qunlun Shen
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy
#from sklearn.metrics import adjusted_mutual_info_score
from autogluon.tabular import TabularPredictor, TabularDataset

def precluster_scRNA_scST(adata, resolution=0.4):
	'''
	Processes single-cell RNA-seq/scST data and get the pre-clustering results.

	 Input:
        adata (anndata.AnnData): An annotated data matrix where rows correspond to cells and columns to genes.
                                 This should be a raw count matrix typically obtained from scRNA-seq experiments.
        resolution (float, optional): A parameter for the Leiden algorithm that affects the granularity of the 
                                      resulting clusters. Higher values lead to more clusters. Default is 0.4.
 	Output:
        adata (anndata.AnnData): The input AnnData object updated with preprocessing and clustering results.
	'''
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	sc.pp.highly_variable_genes(adata, n_top_genes=3000)
	adata = adata[:, adata.var.highly_variable]
	sc.pp.scale(adata, max_value=10)
	sc.pp.pca(adata, n_comps=50)
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
	sc.tl.leiden(adata, resolution=resolution)
	adata.obs['preprocessed'] = True

	return adata

def precluster_scATAC(adata, resolution=0.7):
	'''
	Processes single-cell ATAC data and get the pre-clustering results.

	 Input:
        adata (anndata.AnnData): An annotated data matrix where rows correspond to cells and columns to genomic features.
                                 Typically, this should be a raw count matrix of ATAC-seq peaks.
        resolution (float, optional): A parameter for the Leiden algorithm that affects the granularity of the 
                                      resulting clusters. Higher values lead to more clusters. Default is 0.4.
 	Output:
        adata (anndata.AnnData): The input AnnData object updated with preprocessing and clustering results.
	'''
	import snapatac2 as snap
	snap.pp.select_features(adata, n_features=250000)
	snap.tl.spectral(adata)
	snap.pp.knn(adata)
	snap.tl.leiden(adata, resolution=resolution)

	return adata

