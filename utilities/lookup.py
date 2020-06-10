# Adapted from code written by Daniel Clarke

import pandas as pd
from tqdm import tqdm
import numpy as np
import datetime


def get_dictionary(path):
    '''
    Returns two dictionaries, the first a mapping from symbols to approved gene
    symbols (synonyms), the second a mapping from approved symbols to Entrez
    Gene IDs from an NCBI gene info source designated by path.
    '''
    ncbi = pd.read_csv(path, sep='\t', usecols=[
                       'Symbol', 'Synonyms', 'GeneID'])

    def split_list(v): return v.split('|') if type(v) == str else []
    ncbi['Synonyms'] = ncbi['Synonyms'].apply(split_list)

    # Map existing entities to NCBI Genes
    symbol_lookup = {}
    geneid_lookup = {}
    for i in ncbi.index:
        approved_sym = ncbi.loc[i, 'Symbol']
        v = ncbi.loc[i, 'GeneID']
        geneid_lookup[approved_sym] = v if v != '-' else np.nan

        for sym in [approved_sym] + ncbi.loc[i, 'Synonyms']:
            if not (sym == '-'):
                symbol_lookup[sym] = approved_sym

    return symbol_lookup, geneid_lookup


def get_lookups():
    '''
    Returns two dictionaries, the first a mapping from symbols to approved gene
    symbols (synonyms), the second a mapping from approved symbols to Entrez
    Gene IDs.

    Sources are currently human, mouse, and chicken, in priority order.
    '''
    sources = [
        'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Non-mammalian_vertebrates/Gallus_gallus.gene_info.gz',
        'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz',
        'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
    ]

    symbol_lookup = {}
    geneid_lookup = {}

    for source in tqdm(sources, desc='Gathering sources'):
        sym, gene = get_dictionary(source)
        symbol_lookup.update(sym)
        geneid_lookup.update(gene)

    return symbol_lookup, geneid_lookup


def save_lookup(symbol_lookup, geneid_lookup):
    '''
    Save the lookups as pandas DataFrames. They are saved as:
        symbol_lookup_<year>_<month>.csv
        geneid_lookup_<year>_<month>.csv
    '''
    date = str(datetime.date.today())[0:7].replace('-', '_')
    symbol = pd.DataFrame.from_dict(
        symbol_lookup, orient='index', columns=['Approved Symbol'])
    geneid = pd.DataFrame.from_dict(
        geneid_lookup, orient='index', columns=['Entrez Gene ID'])

    symbol.to_csv('symbol_lookup_{}.csv'.format(date), sep='\t')
    geneid.to_csv('geneid_lookup_{}.csv'.format(date), sep='\t')


def load_lookup(symbol_path, geneid_path):
    '''
    Loads the lookups from custom paths. The files should be tab-separated 
    files. Returns the symbol and geneid lookups as dictionaries.
    '''
    symbol = pd.read_csv(symbol_path, sep='\t', na_filter=False)
    geneid = pd.read_csv(geneid_path, sep='\t', na_filter=False)

    symbol_lookup = dict(zip(symbol.iloc[:, 0], symbol.iloc[:, 1]))
    geneid_lookup = dict(zip(geneid.iloc[:, 0], geneid.iloc[:, 1]))

    return symbol_lookup, geneid_lookup
