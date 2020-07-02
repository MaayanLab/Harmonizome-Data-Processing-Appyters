# %% [markdown]
# # The Cancer Gene Atlas (TCGA)
# %% [markdown]
# Created by: Charles Dai <br>
# Credit to: Moshe Silverstein
# 
# Data Source: https://gdc.cancer.gov/
# %%
# appyter init
from appyter import magic
magic.init(lambda _=globals: _())
# %%
import sys
import os
from datetime import date
import gzip
import io
import requests
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline
from tqdm import tqdm

import harmonizome.utility_functions as uf
import harmonizome.lookup as lookup
# %%
# from clustergrammer_widget import *
# net = Network(clustergrammer_widget)
# %%
%load_ext autoreload
%autoreload 2
# %% [markdown]
# ### Notebook Information
# %%
print('This notebook was run on:', date.today(), '\nPython version:', sys.version)
# %% [markdown]
# # Initialization
# %%
%%appyter hide_code
{% do SectionField(
    name='data',
    title='Load Data',
    subtitle='Upload Files from the Comparative Toxicogenomics Database',
) %}
project_name = 'Mixed germ cell tumor'
# %%
base_url = 'https://api.gdc.cancer.gov/'
files_endpt = base_url + 'files/'
cases_endpt = base_url + 'cases/'
# %% [markdown]
# ### Load Mapping Dictionaries
# %%
symbol_lookup, geneid_lookup = lookup.get_lookups()
# %% [markdown]
# ### Output Path
# %%
output_name = 'tcga'

path = 'Output/TCGA'
if not os.path.exists(path):
    os.makedirs(path)
# %% [markdown]
# # Load Data
# %% [markdown]
# ## File Name and IDs
# %%
fields = [
    'cases.case_id'
]

filters = {
    'op': 'and',
    'content': [{
        'op': 'in',
        'content': {
            'field': 'experimental_strategy',
            'value': ['RNA-Seq'],
        }
    }, 
    {
        'op': 'in',
        'content': {
            'field': 'access',
            'value': ['open'],
        }
    },
    {
        'op': 'in',
        'content': {
            'field': 'file_name',
            'value': ['*htseq.counts.gz'],
        }
    },
    {
        'op': 'in',
        'content': {
            'field': 'cases.project.name',
            'value': [project_name],
        }
    }
    ],
}
# %%
params = {
    'fields': ','.join(fields),
    'filters': json.dumps(filters),
    'size': 100000,
    'facets': 'cases.case_id'
}
response = requests.get('https://api.gdc.cancer.gov/files', params=params)
data = response.json()['data']['hits']
# %%
files = pd.DataFrame([(f['id'], f['cases'][0]['case_id']) for f in data], columns=['file_id', 'case_id']).set_index('file_id')
print(files.shape)
files.head()
# %%
matrix = pd.DataFrame()

for file_id in tqdm(files.index, unit='samples'):
    response = requests.get(data_endpt + file_id, headers = {"Content-Type": "application/json"})
    string_data = io.StringIO(str(gzip.decompress(response.content), 'utf-8'))
    matrix = pd.concat([matrix, pd.read_csv(string_data, sep='\t', header=None, names=['ENSMBL ID', files.loc[file_id, 'case_id']], index_col=0)], axis=1)
matrix.head()
# %%
matrix.index = matrix.index.map(lambda x: x.split('.')[0])
matrix.head()
# %%
cases_fields = requests.get(cases_endpt + '_mapping').json()['fields']
keyfields = [field for field in cases_fields if 
    any(word in field for word in ['demographic', 'diagnoses']) and 'treatment' not in field]
# %%
sample_meta = pd.DataFrame()
for case_id in tqdm(files['case_id'].drop_duplicates(), unit='cases'):
    response = requests.get(cases_endpt + case_id, params={'fields': ','.join(keyfields)})
    data = response.json()['data']
    sample = pd.DataFrame([{'case_id': case_id, **data['demographic'], **data['diagnoses'][0]}])
    sample_meta = pd.concat([sample_meta, sample])
sample_meta = sample_meta.set_index('case_id')
sample_meta.head()
# %% [markdown]
# # Pre-process Data
# %% [markdown]
# ## Create Matrix of Inference Scores
# %%
# Remove duplicates by averaging first
matrix = matrix.set_index(['# GeneSymbol', 'DiseaseName'])
matrix = matrix.groupby(level=[0,1]).mean().unstack()
matrix.head()
# %%
matrix.columns = matrix.columns.droplevel(0)
matrix.columns.name = 'Disease Name'
matrix.index.name = 'Gene Symbol'
matrix.head()
# %%
matrix.shape
# %% [markdown]
# ## Save Unfiltered Matrix to file
# %%
uf.save_data(matrix, path, output_name + '_matrix_unfiltered',
            compression='gzip', dtype=np.float32)
# %% [markdown]
# # Filter Data
# %% [markdown]
# ## Map Gene Symbols to Up-to-date Approved Gene Symbols
# %%
matrix = uf.map_symbols(matrix, symbol_lookup)
matrix.shape
# %% [markdown]
# ## Merge Duplicate Genes By Rows and Duplicate Columns
# %%
matrix = uf.merge(matrix, 'row')
matrix = uf.merge(matrix, 'column')
matrix.shape
# %% [markdown]
# ## Remove Genes that are More Than 95% Missing or Zero Inference
# %%
matrix = matrix.replace(0, np.nan).dropna(
    thresh=0.05 * matrix.shape[1], axis=0).replace(np.nan, 0)
matrix.head()
# %%
matrix.shape
# %% [markdown]
# ## Normalize Matrix (Quantile Normalize the Matrix by Column)
# %%
matrix = uf.quantile_normalize(matrix)
matrix.head()
# %% [markdown]
# ## Normalize Matrix (Z-Score the Rows)
# %%
matrix = uf.zscore(matrix)
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
uf.save_data(matrix, path, output_name + '_matrix_filtered', 
            ext='tsv', compression='gzip')
# %% [markdown]
# # Analyze Data
# %% [markdown]
# ## Create Gene List
# %%
gene_list = uf.gene_list(matrix, geneid_lookup)
gene_list.head()
# %%
gene_list.shape
# %%
uf.save_data(gene_list, path, output_name + '_gene_list',
            ext='tsv', compression='gzip', index=False)
# %% [markdown]
# ## Create Attribute List
# %%
attribute_list = uf.attribute_list(matrix)
attribute_list.head()
# %%
attribute_list.shape
# %%
uf.save_data(attribute_list, path, output_name + '_attribute_list',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Create matrix of Standardized values (values between -1, and 1)
# %%
standard_matrix = uf.standardized_matrix(matrix)
standard_matrix.head()
# %%
uf.save_data(standard_matrix, path, output_name + '_standard_matrix',
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
ternary_matrix = uf.ternary_matrix(standard_matrix)
ternary_matrix.head()
# %%
uf.save_data(ternary_matrix, path, output_name + '_ternary_matrix',
            ext='tsv', compression='gzip')
# %% [markdown]
# ## Create Gene and Attribute Set Libraries
# %%
uf.save_setlib(ternary_matrix, 'gene', 'up', path, output_name + '_gene_up_set')
# %%
uf.save_setlib(ternary_matrix, 'gene', 'down', path, output_name + '_gene_down_set')
# %%
uf.save_setlib(ternary_matrix, 'attribute', 'up', path, 
                           output_name + '_attribute_up_set')
# %%
uf.save_setlib(ternary_matrix, 'attribute', 'down', path, 
                             output_name + '_attribute_down_set')
# %% [markdown]
# ## Create Attribute Similarity Matrix
# %%
attribute_similarity_matrix = uf.similarity_matrix(standard_matrix.T, 'cosine')
attribute_similarity_matrix.head()
# %%
uf.save_data(attribute_similarity_matrix, path,
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
gene_similarity_matrix = uf.similarity_matrix(standard_matrix, 'cosine')
gene_similarity_matrix.head()
# %%
uf.save_data(gene_similarity_matrix, path, 
            output_name + '_gene_similarity_matrix',
            compression='npz', symmetric=True, dtype=np.float32)
# %% [markdown]
# ## Create Gene-Attribute Edge List
# %%
edge_list = uf.edge_list(standard_matrix)
uf.save_data(edge_list, path, output_name + '_edge_list', 
        ext='tsv', compression='gzip')
# %% [markdown]
# # Create Downloadable Save File
# %%
uf.archive(path)
# %% [markdown]
# ### Link to download output files: [click here](./output_archive.zip)
