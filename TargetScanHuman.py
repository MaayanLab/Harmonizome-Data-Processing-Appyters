# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %% [markdown]
# # TargetScanHuman
# %% [markdown]
# Author: Moshe Silverstein <br/>
# Date: 11-17 <br/>
# Data Source: http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi
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
import zipfile

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
output_name = 'tsh'

path = 'Output/TSH'
if not os.path.exists(path):
    os.makedirs(path)
# %%
%%appyter hide_code
{% do SectionField(
    name='data',
    title='Load Data',
    subtitle='Upload Files from the Target Scan Human Prediction Data Sets'
) %}
# %% [markdown]
# # Load Data
# %%
%%appyter code_exec

with zipfile.ZipFile({{FileField(
    constraint='.*\.zip$',
    name='conserved_site_scores', 
    label='Conserved Site Scores', 
    default='Input/TSH/Conserved_Site_Context_Scores.txt.zip',
    section='data')
}}) as zipf:
    with zipf.open(zipf.namelist()[0]) as f:
        conserved = pd.read_csv(f, sep='\t', 
                                usecols=['Gene Symbol', 'miRNA', 'Gene Tax ID'])

# %%
conserved.head()
# %%
conserved.shape
# # Load Data
# %%
%%appyter code_exec

with zipfile.ZipFile({{FileField(
    constraint='.*\.zip$',
    name='nonconserved_site_scores', 
    label='Conserved Site Scores', 
    default='Input/TSH/Nonconserved_Site_Context_Scores.txt.zip',
    section='data')
}}) as zipf:
    with zipf.open(zipf.namelist()[0]) as f:
        nonconserved = pd.read_csv(f, sep='\t',
                                usecols=['Gene Symbol', 'miRNA', 'Gene Tax ID'])

# %%
nonconserved.head()
# %%
nonconserved.shape
# %% [markdown]
# # Pre-process Data
# %% [markdown]
# ## Combine Datasets
# %%
df = pd.concat([conserved, nonconserved])
df.shape
# %% [markdown]
# ## Get Relevant Data
# %%
# Get only human and mouse species
df = df[np.logical_or(df['Gene Tax ID'] == 9606, df['Gene Tax ID'] == 10090)]
df = df.drop('Gene Tax ID', axis=1).set_index('Gene Symbol')
df.head()
# %%
df.shape
# %% [markdown]
# # Filter Data
# %% [markdown]
# ## Map Gene Symbols to Up-to-date Approved Gene Symbols
# %%
df = uf.mapgenesymbols(df, symbol_lookup)
df.shape
# %% [markdown]
# # Analyze Data
# %% [markdown]
# ## Create Binary Matrix
# %%
binary_matrix = uf.createBinaryMatrix(df)
binary_matrix.head()
# %%
binary_matrix.shape
# %%
uf.saveData(binary_matrix, path, output_name + '_binary_matrix', 
            compression='npz', dtype=np.uint8)
# %% [markdown]
# ## Create Gene List
# %%
gene_list = uf.createGeneList(binary_matrix, geneid_lookup)
gene_list.head()
# %%
gene_list.shape
# %%
uf.saveData(gene_list, path, output_name + '_gene_list',
            ext='tsv', compression='gzip', index=False)
# %% [markdown]
# ## Create Attribute List
# %%
attribute_list = uf.createAttributeList(binary_matrix)
attribute_list.head()
# %%
attribute_list.shape
# %%
uf.saveData(attribute_list, path, output_name + '_attribute_list',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Create Gene and Attribute Set Libraries
# %%
uf.createUpGeneSetLib(binary_matrix, path, output_name + '_gene_up_set')
# %%
uf.createUpAttributeSetLib(binary_matrix, path, 
                           output_name + '_attribute_up_set')
# %% [markdown]
# ## Create Attribute Similarity Matrix
# %%
attribute_similarity_matrix = uf.createSimilarityMatrix(binary_matrix.T, 'jaccard', sparse=True)
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
gene_similarity_matrix = uf.createSimilarityMatrix(binary_matrix, 'jaccard', sparse=True)
gene_similarity_matrix.head()
# %%
uf.saveData(gene_similarity_matrix, path, 
            output_name + '_gene_similarity_matrix',
            compression='npz', symmetric=True, dtype=np.float32)
# %% [markdown]
# ## Create Gene-Attribute Edge List
# %%
uf.createGeneAttributeEdgeList(binary_matrix, attribute_list, gene_list, 
                               path, output_name + '_gene_attribute_edge_list')
# %% [markdown]
# # Create Downloadable Save File
# %%
uf.createArchive(path)
# %% [markdown]
# ### Link to download output files: [click here](./output_archive.zip)
