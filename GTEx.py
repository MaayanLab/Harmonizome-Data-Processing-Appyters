# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # GTEx Tissue
# %% [markdown]
# Author: Moshe Silverstein <br/>
# Date: 12-2017 <br/>
# Data Source: https://www.gtexportal.org/home/datasets
#
# Reviewer: Charles Dai <br>
# Updated: 6-20
# %%
# appyter init
from appyter import magic
magic.init(lambda _=globals: _())
# %%
import sys
import os
import gzip

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

import utility_functions as uf
import lookup
# %%
# from clustergrammer_widget import *
# net = Network(clustergrammer_widget)
# %%
%load_ext autoreload
%autoreload 2
# %% [markdown]
# ### Python Version
# %%
sys.version
# %% [markdown]
# # Initialization
# %% [markdown]
# ### Load Mapping Dictionaries
# %%
symbol_lookup, geneid_lookup = lookup.get_lookups()
# %% [markdown]
# ### Output Path
# %%
output_name = 'gtex'

path = 'Output/GTEx'
if not os.path.exists(path):
    os.makedirs(path)
# %%
%%appyter hide_code
{% do SectionField(
    name='data',
    title='Load Data',
    subtitle='Upload Files from the GTEx Portal',
) %}
# %% [markdown]
# # Load Data
# %%
%%appyter code_exec

matrix_zip = {{FileField(
    constraint='.*\.gz$',
    name='expression_matrix', 
    label='RNA-Seq Gene TPMs (gz)', 
    default='Input/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz',
    section='data')
}}
with gzip.open(matrix_zip, 'r') as matrix_file:
    matrix = pd.read_csv(matrix_file, sep='\t', skiprows=2)
# %%
matrix.head()
# %%
matrix.shape
# %% [markdown]
# ## Load Sample Metadata
# %%
%%appyter code_exec

sample_meta = pd.read_csv({{FileField(
    constraint='.*\.txt$',
    name='sample_metadata', 
    label='Sample Metadata', 
    default='Input/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
    section='data')
}}, sep='\t', index_col=0, usecols='SMTSD')
# %%
sample_meta.head()
# %%
sample_meta.shape
'''
# %% [markdown]
# ## Load Subject Metadata
# %%
%%appyter code_exec

subject_meta = pd.read_csv({{FileField(
    constraint='.*\.txt$',
    name='subject_metadata', 
    label='Subject Metadata', 
    default='Input/GTEx/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt',
    section='data')
}})
# %%
subject_meta.head()
# %%
subject_meta.shape
'''
# %% [markdown]
# # Pre-process Data
# %% [markdown]
# ## Map Sample ID to Attribute
# %%
matrix.columns = sample_meta.reindex(matrix.columns).reset_index(drop=True)
matrix.head()
# %% [markdown]
# ## Save Unfiltered Matrix to file
# %%
uf.saveData(matrix, path, output_name + '_matrix_unfiltered',
            compression='gzip', dtype=np.float32)
# %% [markdown]
# # Filter Data
# %% [markdown]
# ## Remove Data that is More Than 95% Missing and Impute Missing Data
# %%
matrix = uf.removeAndImpute(matrix)
matrix.head()
# %%
matrix.shape
# %% [markdown]
# ## Map Gene Symbols to Up-to-date Approved Gene Symbols
# %%
matrix = uf.mapgenesymbols(matrix, symbol_lookup)
matrix.shape
# %% [markdown]
# ## Merge Duplicate Genes By Rows and Duplicate Columns
# %%
matrix = uf.merge(matrix, 'row', 'mean')
matrix = uf.merge(matrix, 'column', 'mean')
matrix.shape
# %% [markdown]
# ## Log2 Transform
# %%
matrix = uf.log2(matrix)
matrix.head()
# %% [markdown]
# ## Normalize Matrix (Quantile Normalize the Matrix by Column)
# %%
matrix = uf.quantileNormalize(matrix)
matrix.head()
# %% [markdown]
# ## Normalize Matrix (Z-Score the Rows)
# %%
matrix = uf.zscore(matrix, 'row')
matrix.head()
# %% [markdown]
# ## Histogram of First Sample
# %%
matrix.iloc[:, 0].hist(bins=100)
# %% [markdown]
# ## Histogram of First Gene
# %%
matrix.iloc[0, :].hist(bins=100)
# %% [markdown]
# ## Save Filtered Matrix
# %%
uf.saveData(matrix, path, output_name + '_matrix_filtered', 
            ext='tsv', compression='gzip')
# %% [markdown]
# # Analyze Data
# %% [markdown]
# ## Create Gene List
# %%
gene_list = uf.createGeneList(matrix, geneid_lookup)
gene_list.head()
# %%
gene_list.shape
# %%
uf.saveData(gene_list, path, output_name + '_gene_list',
            ext='tsv', compression='gzip', index=False)
# %% [markdown]
# ## Create Attribute List
# %%
attribute_list = uf.createAttributeList(matrix)
attribute_list.head()
# %%
attribute_list.shape
# %%
uf.saveData(attribute_list, path, output_name + '_attribute_list',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Create matrix of Standardized values (values between -1, and 1)
# %%
standard_matrix = uf.createStandardizedMatrix(matrix)
standard_matrix.head()
# %%
uf.saveData(standard_matrix, path, output_name + '_standard_matrix',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Plot of A Single Celltype, Normalized Value vs. Standardized Value
# %%
plt.plot(matrix[matrix.columns[0]],
         standard_matrix[standard_matrix.columns[0]], 'bo')
plt.xlabel('Normalized Values')
plt.ylabel('Standardized Values')
plt.title(standard_matrix.columns[0])
plt.grid(True)
# %% [markdown]
# ## Create Ternary Matrix
# %%
ternary_matrix = uf.createTernaryMatrix(standard_matrix)
ternary_matrix.head()
# %%
uf.saveData(ternary_matrix, path, output_name + '_ternary_matrix',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Create Gene and Attribute Set Libraries
# %%
uf.createUpGeneSetLib(ternary_matrix, path, output_name + '_gene_up_set')
# %%
uf.createDownGeneSetLib(ternary_matrix, path, output_name + '_gene_down_set')
# %%
uf.createUpAttributeSetLib(ternary_matrix, path, 
                           output_name + '_attribute_up_set')
# %%
uf.createDownAttributeSetLib(ternary_matrix, path, 
                             output_name + '_attribute_down_set')
# %% [markdown]
# ## Create Attribute Similarity Matrix
# %%
attribute_similarity_matrix = uf.createSimilarityMatrix(matrix.T, 'cosine')
attribute_similarity_matrix.head()
# %%
uf.saveData(attribute_similarity_matrix, path,
            output_name + '_attribute_similarity_matrix', 
            compression='npz', symmetric=True, dtype=np.float32)
# %%
# net.load_df(attribute_similarity_matrix.iloc[:,:].copy())
# net.filter_N_top('row', rank_type='sum', N_top=300)
# net.cluster()
# net.widget()
# %% [markdown]
# ## Create Gene Similarity Matrix
# %%
gene_similarity_matrix = uf.createSimilarityMatrix(matrix, 'cosine')
gene_similarity_matrix.head()
# %%
uf.saveData(gene_similarity_matrix, path, 
            output_name + '_gene_similarity_matrix',
            compression='npz', symmetric=True, dtype=np.float32)
# %% [markdown]
# ## Create Gene-Attribute Edge List
# %%
uf.createGeneAttributeEdgeList(standard_matrix, attribute_list, gene_list, 
                               path, output_name + '_gene_attribute_edge_list')
# %% [markdown]
# # Create Downloadable Save File
# %%
uf.createArchive(path)
# %% [markdown]
# ### Link to download output files: [click here](./output_archive.zip)