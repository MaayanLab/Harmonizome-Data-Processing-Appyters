# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# # Allen Brain Atlas Developing Human Brain Sample (microarray)
# %% [markdown]
# Author: Moshe Silverstein <br/>
# Date: 3-18 <br/>
# Data Source: http://www.brainspan.org/static/download.html
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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

import harmonizome.utility_functions as uf
import harmonizome.lookup as lookup
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
output_name = 'aba_dev_mic_samp'

path = 'Output/ABA-DHB-MA-S'
if not os.path.exists(path):
    os.makedirs(path)
# %%
%%appyter hide_code
{% do SectionField(
    name='data',
    title='Load Data',
    subtitle='Upload Files from the Developing Human Brain Sample (microarray) Data Set. Link: http://www.brainspan.org/static/download.html',
) %}
# %% [markdown]
# # Load Data
# %%
%%appyter code_exec

matrix = pd.read_csv({{FileField(
    constraint='.*\.csv$',
    name='expression_matrix', 
    label='Expression Matrix', 
    default='Input/ABA-DHB-MA-S/gene_array_matrix_csv/expression_matrix.csv',
    section='data')
}}, header=None, index_col=0)
# %%
matrix.head()
# %%
matrix.shape
# %% [markdown]
# ## Load Sample Metadata
# %%
%%appyter code_exec

sample_meta = pd.read_csv({{FileField(
    constraint='.*\.csv$',
    name='columns_metadata', 
    label='Sample Metadata', 
    default='Input/ABA-DHB-MA-S/gene_array_matrix_csv/columns_metadata.csv',
    section='data')
}}, index_col=0)
# %%
sample_meta.head()
# %%
sample_meta.shape
# %% [markdown]
# ## Load Gene Metadata
# %%
%%appyter code_exec

gene_meta = pd.read_csv({{FileField(
    constraint='.*\.csv$',
    name='rows_metadata', 
    label='Gene Metadata', 
    default='Input/ABA-DHB-MA-S/gene_array_matrix_csv/rows_metadata.csv',
    section='data')
}}, index_col=0)
# %%
gene_meta.head()
# %%
gene_meta.shape
# %% [markdown]
# # Pre-process Data
# %% [markdown]
# ## Map Sample Metadata to Sample ID
# %%
matrix.columns = sample_meta.apply(lambda row: '_'.join(
    [row['structure_name'], str(row['age']), row['gender'], str(row['donor_id'])]), axis=1)
matrix.head()
# %% [markdown]
# ## Map Gene to Row
# %%
matrix.index = gene_meta['gene_symbol']
matrix.index.name = None
matrix.head()
# %%
matrix.shape
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
