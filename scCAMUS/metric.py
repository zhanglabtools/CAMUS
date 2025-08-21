# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 20:52:22 2024

@author: Qunlun Shen
"""

from sklearn.metrics import adjusted_mutual_info_score
from autogluon.tabular import TabularPredictor, TabularDataset
from scipy.stats import entropy
import numpy as np
import pandas as pd
import scanpy as sc
from .preprocess import precluster_scRNA_scST, precluster_scATAC

def CAMUS_score(annotations, pre_clusters):
	'''
	Calculate the CAMUS score to evaluate the coaccordance between the predicted annotation results and the preclustered results.

	Parameters:
		annotations (list or array-like): The annotation results for the dataset, should have the same length as `pre_clusters`.
		pre_clusters (list or array-like): The pre-cluster labels, should have the same length as `annotations`.

	Returns:
		float: The CAMUS Score representing the quality of the clustering.

	Raises:
		ValueError: If the input arrays `annotations` and `pre_clusters` do not have the same length.
	'''
	if len(annotations) != len(pre_clusters):
		raise ValueError("The length of 'annotations' and 'pre_clusters' must be the same.")

	score = adjusted_mutual_info_score(annotations, pre_clusters)

	return score


def CAMUS_prioritize(adata, key_class=['annot']):
	'''
	Processes single-cell data to get the pre-clustering results and 
	prioritizes the annotation results.

	 Input:
        adata (anndata.AnnData): The annotated data matrix where rows correspond to cells and columns to features.
        key_class (list of str, optional): List of keys in `adata.obs` to use for annotation comparison.
                                           Defaults to ['annot'].
        resolution (float, optional): Clustering resolution for the Leiden algorithm, affecting cluster granularity.
                                      Defaults to 0.4.
        data_type (str, optional): Specifies the type of data processing to apply, supporting 'scRNA-seq', 'scST',
                                   or 'scATAC-seq'. Defaults to 'scRNA-seq'.

    Output:
        dict: A dictionary where keys are the annotation keys from `key_class` and values are their corresponding
              CAMUS scores, sorted from highest to lowest based on CAMUS scores.
	'''
	#if data_type in ['scRNA-seq', 'scST']:
		#adata = precluster_scRNA_scST(adata, resolution=resolution)
	#if data_type in ['scATAC-seq']:
		#adata = precluster_scATAC(adata, resolution=resolution)

	#adata.obs['preprocessed'] = True
	dict_prioritize = {}
	pre_clusters = adata.obs['leiden'].tolist()
	for annot_key in key_class:
		annotation = adata.obs[annot_key].tolist()
		score = CAMUS_score(annotation, pre_clusters)
		dict_prioritize[annot_key] = score

	sorted_dict = sorted(dict_prioritize.items(), key=lambda item: item[1], reverse=True)
	sorted_dict = dict(sorted_dict)

	return sorted_dict

def CAMUS_estimate_acc(adata, annot_key='annot', data_type='scRNA-seq', model_path=None):
	"""
	Estimates the accuracy of annotations in single-cell data using a trained machine learning model.

	Inputs:
	adata: AnnData object containing single-cell data. This object must have a 'preprocessed' column in `adata.obs` indicating prior preprocessing.
	key_class: str, default 'annot'. Name of the column in `adata.obs` containing the annotations for cells.
	data_type: str, default 'scRNA-seq'. Type of single-cell data. Valid options: 'scRNA-seq', 'scST', 'scATAC-seq'.

	Output:
	estimated_accuracy: The estimated accuracy of the annotations.

    """
	if 'preprocessed' not in adata.obs.columns:
		raise ValueError('Please run CAMUS_prioritize first!')
	if model_path==None:
		raise ValueError('Please provide the path of the trained model!')

	predictor = TabularPredictor.load(path=model_path)
	for resolution in [0.2, 0.4, 0.6, 0.8]:
		if data_type in ['scRNA-seq', 'scST']:
			sc.tl.leiden(adata, resolution=resolution, key_added='leiden_'+str(resolution))
		if data_type in ['scATAC-seq']:
			import snapatac2 as snap
			snap.tl.leiden(adata, resolution=resolution, key_added='leiden_'+str(resolution))
	
	X_test = [adata.shape[0]]
	annotations = adata.obs[annot_key].tolist()
	for resolution in [0.2, 0.4, 0.6, 0.8]:
		pre_clusters = adata.obs['leiden_'+str(resolution)].tolist()
		cluster_num = len(set(pre_clusters))
		cluster_counts = np.bincount(pre_clusters, minlength=len(set(pre_clusters)))
		#cluster_sizes = np.mean(cluster_counts)
		probabilities = cluster_counts / cluster_counts.sum()
		cluster_entropy = entropy(probabilities)
		score = CAMUS_score(annotations, pre_clusters)
		X_test.append(score)
		X_test.append(cluster_num)
		X_test.append(cluster_entropy)
	estimated_accuracy = predictor.predict(pd.DataFrame(X_test, index=['0', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13']).T)
	return estimated_accuracy
