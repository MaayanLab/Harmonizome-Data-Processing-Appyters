# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %% [markdown]
# # Guide to Pharmacology (Chemical Ligands of Receptors)
# %% [markdown]
# Author: Moshe Silverstein <br/>
# Date: 8-17 <br/>
# Data Source: http://www.guidetopharmacology.org/download.jsp
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
# ### Choose Interaction
%%appyter code_eval
{% set ligand = ChoiceField(
    name='ligands',
    label='Ligands of receptors',
    choices={
        'Protein Ligands': "'protein'",
        'Chemical Ligands': "'chemical'"
    },
    default='Protein Ligands',
    section='data'
) %}
# %% [markdown]
# ### Load Mapping Dictionaries
# %%
symbol_lookup, geneid_lookup = lookup.get_lookups()
# %% [markdown]
# ### Output Path
# %%
%%appyter code_exec

output_name = 'guide_to_pharm_' + {{ligand}}

path = 'Output/Guide-to-Pharm-' + {{ligand}}.upper()
if not os.path.exists(path):
    os.makedirs(path)
# %%
%%appyter hide_code
{% do SectionField(
    name='data',
    title='Load Data',
    subtitle='Upload Files from the Guide to Pharmacology Interactions Data Set',
) %}
# %% [markdown]
# # Load Data
# %%
%%appyter code_exec

headers = ['target_species', 'ligand_gene_symbol', 'target_gene_symbol', 
            'ligand']
df = pd.read_csv({{FileField(
    constraint='.*\.csv$',
    name='interactions', 
    label='Interaction Data', 
    default='Input/Guide-to-Pharm/interactions.csv',
    section='data')
}}, usecols=headers)
# %%
df.head()
# %%
df.shape
# %% [markdown]
# # Pre-process Data
# %% [markdown]
# ## Get Relevant Data
# %%
df = df[df['target_species'].isin(['Human', 'Mouse', 'Rat'])]]
# %%
%%appyter code_exec

if {{ligand}} == 'protein':
    df = df[['target_gene_symbol', 'ligand_gene_symbol']]
elif {{ligand}} == 'chemical':
    df = df[['target_gene_symbol', 'ligand']]
df = df.dropna()
df['target_gene_symbol'] = df['target_gene_symbol'].str.split('|')
df = df.explode('target_gene_symbol').set_index('target_gene_symbol')
df.head()
# %% [markdown]
# # Filter Data
# %% [markdown]
# ## Map Gene Symbols to Up-to-date Approved Gene Symbols
# %%
%%appyter code_exec
if {{ligand}} == 'protein':
    df = df.reset_index().set_index('ligand_gene_symbol')
    df = uf.mapgenesymbols(df, symbol_lookup)
    df = df.reset_index().set_idnex('target_gene_symbol')
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
attribute_similarity_matrix = uf.createSimilarityMatrix(binary_matrix.T, 'jaccard', np.bool)
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
gene_similarity_matrix = uf.createSimilarityMatrix(binary_matrix, 'jaccard', np.bool)
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
