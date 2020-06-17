# Adapted from code created by Moshe Silverstein

import datetime
import os
import zipfile

import numpy as np
import pandas as pd
import scipy.spatial.distance as dist
import scipy.sparse as sp
from statsmodels.distributions.empirical_distribution import ECDF

from tqdm import tqdm


def removeAndImpute(inputDF):
    '''
    Removes rows and columns that have more than 95% of their data missing,
    i.e. set to 0. Replacing any missing data leftover after removal with
    the means of the rows.
    '''
    outputDF = inputDF.replace(0.0, np.nan)
    outputDF = outputDF.dropna(thresh=0.05 * inputDF.shape[0], axis=1)
    outputDF = outputDF.dropna(thresh=0.05 * inputDF.shape[1], axis=0)

    outputDF = outputDF.T.fillna(
        outputDF.mean(axis=1, numeric_only=True), axis=0).T
    return outputDF


def merge(inputDF, axis, method):
    '''
    Merges duplicate rows or columns, depending on the axis specified. The
    final values of the merged rows or columns is determined by the method.
    '''
    axis = axis.lower()
    if method == 'mean':
        if axis == 'column':
            outputDF = inputDF.groupby(inputDF.columns, axis=1).mean()
        elif axis == 'row':
            outputDF = inputDF.groupby(level=0, axis=0).mean()
        return outputDF


def quantileNormalize(inputDF):
    '''
    Performs quantile normalization on the input DataFrame.
    '''
    # from ayhan on StackOverflow
    rank_mean = inputDF.stack().groupby(
        inputDF.rank(method='first').stack().astype(int)).mean()
    outputDF = inputDF.rank(method='min').stack().astype(int).map(
        rank_mean).unstack()
    return outputDF


def zscore(inputDF, axis, epsilon=0):
    '''
    Calculates the modified z-score of inputDF according to the specified axis.

    Parameters:
        axis - the axis on which to calculate the z-scores. Either 'row' or 'column'
        epsilon - small adjustment in the case of divide by 0 errors.
    '''
    np.seterr(divide='ignore', invalid='ignore')
    if axis == 'row':
        a = inputDF.to_numpy()
        median_y = np.median(a, axis=1)[:, np.newaxis]
        abs_dev = np.abs(a - median_y)
        median_dev = np.median(abs_dev, axis=1)
        mean_dev = np.mean(abs_dev, axis=1)
        median_abs_dev = np.broadcast_to(median_dev[:, np.newaxis], a.shape)
        mean_abs_dev = np.broadcast_to(mean_dev[:, np.newaxis], a.shape)
        modified_z_scores = np.where(median_abs_dev != 0,
                                     0.6745 * (a - median_y) / median_abs_dev,
                                     (a - median_y) / (1.253314 * mean_abs_dev + epsilon))

        outputDF = pd.DataFrame(data=modified_z_scores, index=inputDF.index,
                                columns=inputDF.columns)
        return outputDF

    elif axis == 'column':
        # left unfactored. Is this used?
        for i, col in enumerate(tqdm(inputDF.columns)):

            mean = inputDF[col].mean()
            std = inputDF[col].std()
            inputDF[col] = inputDF[col].apply(lambda x: ((x-mean)/std))


def log2(inputDF):
    '''
    Returns a dataframe with the log2 values of the input.
    '''
    a = inputDF.to_numpy(copy=True)
    a = np.log2(a + 1)
    outputDF = pd.DataFrame(data=a, index=inputDF.index,
                            columns=inputDF.columns)
    return outputDF


def mapgenesymbols(inputDF, symbol_lookup):
    '''
    Replaces the index of the inputDF, which are gene names, with
    corresponding approved gene symbols according to the given symbol_lookup 
    dictionary. If any gene names are not in the mapping, they are discarded 
    from the DataFrame.
    '''
    tqdm.pandas()
    inputDF = inputDF.reset_index()

    inputDF.iloc[:, 0] = inputDF.iloc[:, 0].progress_map(
        lambda x: symbol_lookup.get(x, np.nan))

    outputDF = inputDF.dropna(subset=[inputDF.columns[0]])
    outputDF = outputDF.set_index(outputDF.columns[0])
    return outputDF


def createTernaryMatrix(inputDF):
    '''
    Returns the input matrix with all significant values, greater than 0.95
    or less than -0.95, mapped to 1 or -1, respectively. All other values
    are mapped to 0.
    '''
    def mapter(x):
        if x >= 0.95:
            return 1
        elif x <= -0.95:
            return -1
        else:
            return 0

    df = inputDF.applymap(mapter)
    return df


def createSetLibHelper(lib, inputDF, path, name, direction, details=None):
    '''
    If lib = 'gene', this creates a file which lists all attributes and the
    genes that are correlated in the direction given with that attribute.
    Only attributes that are correlated with more than 5 and less than or equal
    to 2000 genes are saved.

    If lib = 'attribute', this creates a file which lists all genes and the
    attributes that are correlated in the direction given with that gene.
    The year and month are added at the end of the name. The path the file is
    saved to is thus
        path + name + '_<year>_<month>.gmt'
    '''
    filenameGMT = getFileName(path, name, 'gmt')

    if not (lib == 'gene' or lib == 'attribute'):
        return
    if lib == 'attribute':
        inputDF = inputDF.T

    with open(filenameGMT, 'w') as f:
        arr = inputDF.reset_index(drop=True).to_numpy(dtype=np.int_)
        attributes = inputDF.columns

        if lib == 'gene':
            num_match = (arr == direction).sum(axis=0)
            sufficient_matched = np.logical_and(
                num_match > 5, num_match <= 2000)
            arr = arr[:, sufficient_matched]
            attributes = attributes[sufficient_matched]

        w, h = arr.shape
        for i in tqdm(range(h)):
            print(attributes[i], details[i] if details else 'NA',
                  *inputDF.index[arr[:, i] == direction],
                  sep='\t', end='\n', file=f)


def createUpGeneSetLib(inputDF, path, name, details=None):
    '''
    Create a file which lists all attributes and the genes that are positively
    correlated with that attribute. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    '''
    createSetLibHelper('gene', inputDF, path, name, 1, details)


def createDownGeneSetLib(inputDF, path, name, details=None):
    '''
    Create a file which lists all attributes and the genes that are negatively
    correlated with that attribute. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    '''
    createSetLibHelper('gene', inputDF, path, name, -1, details)


def createUpAttributeSetLib(inputDF, path, name):
    '''
    Create a file which lists all genes and the attributes that are positively
    correlated with that gene. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    '''
    createSetLibHelper('attribute', inputDF, path, name, 1)


def createDownAttributeSetLib(inputDF, path, name):
    '''
    Create a file which lists all genes and the attributes that are negatively
    correlated with that gene. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    '''
    createSetLibHelper('attribute', inputDF, path, name, -1)


def createSimilarityMatrix(inputDF, metric, dtype=None, sparse=False):
    '''
    Creates a similarity matrix between the rows of the inputDF based on
    the metric specified. The resulting matrix has both rows and columns labeled
    by the index of inputDF.
    '''
    if sparse and metric == 'jaccard':
        # from na-o-ys on Github
        sparse = sp.csr_matrix(inputDF.to_numpy(dtype=np.bool).astype(np.int))
        cols_sum = sparse.getnnz(axis=1)
        ab = sparse * sparse.T
        aa = np.repeat(cols_sum, ab.getnnz(axis=1)) # for rows
        bb = cols_sum[ab.indices] # for columns
        similarities = ab.copy()
        similarities.data = similarities.data / (aa + bb - ab.data)
        similarity_matrix = similarities.todense()
        np.fill_diagonal(similarity_matrix, 1)

    else:
        similarity_matrix = dist.pdist(inputDF.to_numpy(dtype=dtype), metric)
        similarity_matrix = dist.squareform(similarity_matrix)
        similarity_matrix = 1 - similarity_matrix

    similarity_df = pd.DataFrame(
        data=similarity_matrix, index=inputDF.index, columns=inputDF.index)
    similarity_df.index.name = None
    similarity_df.columns.name = None
    return similarity_df


def createGeneList(inputDF, geneid_lookup):
    '''
    Creates a list of genes and the corresponding Entrez Gene IDs(supplied by
    the NCBI)

    Note: this differs from the previous function in its behavior with dealing
    with genes that do not have an ID. This function will set the id of the gene
    to -1, whereas the previous script will set them to np.nan.
    '''
    gene_list = inputDF.index
    gene_ids = np.array([geneid_lookup.get(x, -1)
                         if np.isfinite(geneid_lookup.get(x, -1))
                         else -1 for x in tqdm(gene_list)], dtype=np.int_)
    df = pd.DataFrame(list(zip(gene_list, gene_ids)),
                      columns=['GeneSym', 'GeneID'])
    return df


def createAttributeList(inputDF, metaData=None):
    '''
    Creates a list of attributes in the form of a DataFrame, with the attributes
    as the indices. If metaData is specified, it returns appends the attributes
    of inputDF onto the metaData DataFrame.
    '''
    if metaData is not None:
        attribute_list = metaData.set_index(inputDF.columns)
    else:
        attribute_list = pd.DataFrame(index=inputDF.columns)
    attribute_list.index.name = 'Attributes'
    return attribute_list


def createGeneAttributeEdgeList(df, attributelist, genelist, path, name):
    '''
    Creates the gene-attribute edge list from the given input DataFrame,
    attribute and gene lists. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    Also prints the number of cells in df that are statistically
    significant, i.e. > 0.95 confidence.
    Requires:
        attributelist and genelist were generated from running
        createAttributeList and createGeneList on df, respectively.
    '''

    count = np.sum(np.sum(df >= 0.95) + np.sum(df <= -0.95))
    df = df.stack()
    df.index.names = ['Gene', 'Attribute']
    df.name = 'Weight'
    df = df.astype('category')
    #df['Weight'] = df['Weight'].astype('category')

    df.columns = ['Gene', 'Attribute', 'Weight']
    ##saveData(df, path, name, ext='tsv', compression='gzip')

    print('The number of statisticaly relevent gene-attribute associations is: %d' % count)


def createStandardizedMatrix(inputDF):
    '''
    Creates a standardized matrix by using an emperical CDF for each row.
    Each row in the inputDF should represent a single gene.

    Requires:
    Indices of the DataFrame are unique.
    '''
    arr = inputDF.to_numpy(copy=True)

    def process(array):
        ourECDF = ECDF(array)
        array = ourECDF(array)
        mean = np.mean(array)
        array = 2 * (array - mean)
        return array

    for i in tqdm(range(arr.shape[0])):
        arr[i, :] = process(arr[i, :])

    values = arr.flatten()
    ourECDF = ECDF(values)
    ourECDF = ourECDF(values).reshape(arr.shape)

    mean = np.mean(ourECDF)
    ourECDF = 2 * (ourECDF - mean)
    newDF = pd.DataFrame(data=ourECDF, index=inputDF.index,
                         columns=inputDF.columns)
    return newDF


def getFileName(path, name, ext):
    '''
    Returns the file name by taking the path and name, adding the year and month
    and then the extension. The final string returned is thus
        '<path>/<name>_<year>_<month>.ext'
    '''
    date = str(datetime.date.today())[0:7].replace('-', '_')
    filename = ''.join([name, '_', date, '.', ext])
    return os.path.join(path, filename)


def saveData(inputDF, path, name, compression=None, ext='tsv',
             symmetric=False, dtype=None, **kwargs):
    '''
    Save inputDF according to the compression method given. 
    compression can take these values:
        None or 'gmt' - defaults to pandas to_csv() function.
        'gzip' - uses the gzip compression method of the pandas to_csv() function
        'npz' - converts the DataFrame to a numpy array, and saves the array.
                The array is stored as 'axes[0]_axes[1]'. If symmetric is true,
                it is stored as 'axes[0]_axes[1]_symmetric' instead.
    ext is only used if compression is None or 'gzip'. The extension of the file
    will be .ext, or .ext.gz if 'gzip' is specified.
    axes must only be specified if compression is 'npz'. It is a string tuple
    that describes the index and columns inputDF, i.e. (x, y) where x, y = 
    'gene' or 'attribute'.
    symmetric is only used if compression is 'npz', and indicates if inputDF
    is symmetric and can be stored as such. 
    dtype is only used if compression is 'npz', and indicates a dtype that the
    array can be cast to before storing.

    The year and month are added at the end of the name. The path the file is 
    saved to is thus
        path + name + '_<year>_<month>.ext'
    where ext is .ext, .ext.gz, or .npz depending on the compression method.
    '''

    if compression is None:
        name = getFileName(path, name, ext)
        inputDF.to_csv(name, sep='\t', **kwargs)
    elif compression == 'gzip':
        name = getFileName(path, name, ext + '.gz')
        inputDF.to_csv(name, sep='\t', compression='gzip', **kwargs)
    elif compression == 'npz':
        name = getFileName(path, name, 'npz')

        data = inputDF.to_numpy(dtype=dtype)
        index = np.array(inputDF.index)
        columns = np.array(inputDF.columns)

        if symmetric:
            data = np.triu(data)
            np.savez_compressed(name, symmetric=data, index=index)
        else:
            np.savez_compressed(name, nonsymmetric=data, 
                                index=index, columns=columns)

def loadData(filename):
    '''
    Loads a pandas DataFrame stored in a .npz data numpy array format.
    '''
    with np.load(filename, allow_pickle=True) as data_load:
        arrays = data_load.files
        if arrays[0] == 'symmetric':
            data = data_load['symmetric']
            index = data_load['index']
            data = data + data.T - np.diag(data.diagonal())
            df = pd.DataFrame(data=data, index=index, columns=index)
            return df
        elif arrays[0] == 'nonsymmetric':
            data = data_load['nonsymmetric']
            index = data_load['index']
            columns = data_load['columns']
            df = pd.DataFrame(data=data, index=index, columns=columns)
            return df


def createArchive(path):
    with zipfile.ZipFile('output_archive.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(path):
            for f in files:
                zipf.write(os.path.join(root, f))


def createBinaryMatrix(inputDF, ppi=False):
    '''
    Creates an adjacency matrix from inputDF, which is a gene-attribute edge
    list.
    '''
    if ppi:
        # Left unfactored. Is this used?
        genes = list(set(inputDF.iloc[:, 0].unique(
        ).tolist()+inputDF.iloc[:, 1].unique().tolist()))
        matrix = pd.DataFrame(index=genes, columns=genes, data=0)
        for i, gene in enumerate(tqdm(genes)):
            lst = inputDF[inputDF.iloc[:, 0] == gene].iloc[:, 1].tolist()
            lst += inputDF[inputDF.iloc[:, 1] == gene].iloc[:, 0].tolist()
            lst = set(lst)
            lst.discard(gene)
            lst = list(lst)
            matrix.loc[gene, lst] = 1
        return(matrix)
    else:
        matrix = pd.crosstab(inputDF.index, inputDF.iloc[:, 0])
        matrix[matrix > 1] = 1
        matrix = matrix.astype('boolean')
        matrix.index.name = None
        matrix.columns.name = None
        return matrix
