import datetime
import numpy as np
import pandas as pd
import scipy.spatial.distance as dist
from statsmodels.distributions.empirical_distribution import ECDF
from tqdm import tqdm


mappingDF = pd.read_csv('mappingFile_2017.txt',
                        sep='\t', header=None, index_col=0)

getGeneIDs = pd.read_csv('GeneSymbolAndIDS_2017.txt', sep='\t', index_col=0)


def removeAndImpute(inputDF):
    df = inputDF.copy()
    # mean = np.mean(df.values.flatten())
    df.reset_index(inplace=True)
    df.replace(0.0, np.nan, inplace=True)
    df.dropna(thresh=int(0.05*df.shape[0]), axis=1, inplace=True)
    df.dropna(thresh=int(0.05*df.shape[1]), axis=0, inplace=True)
    genes = df.iloc[:, 0].values.tolist()
    df.drop(df.columns[0], axis=1, inplace=True)
    # df.replace(np.nan, mean, inplace=True)
    df = df.T.fillna(df.mean(axis=1)).T
    df.index = genes
    return(df)


def merge(inputDF, axis, method):
    axis = axis.lower()
    if method == 'mean':
        if axis == 'column':
            inputDF = inputDF.groupby(inputDF.columns, axis=1).mean()
        elif axis == 'row':
            inputDF = inputDF.groupby(level=0, axis=0).mean()
        return(inputDF)


def quantileNormalize(inputDF):
    df = inputDF.copy()

    attributes = df.columns.values.tolist()

    df.columns = np.arange(0, len(attributes))

    # compute rank
    dic = {}
    for i, col in enumerate(tqdm(list(df))):

        dic.update({col: sorted(df[col])})

    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis=1).tolist()
    # sort
    for i, col in enumerate(tqdm(list(df))):

        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]

    df.columns = attributes
    return df


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
        inputDF.update(outputDF)

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
    return(outputDF)


def mapgenesymbols(inputDF):
    '''
    Replaces the index of the inputDF, which are gene names, with
    corresponding approved gene symbols from a built-in mapping.
    If any gene names are not in the mapping, they are discarded from the
    DataFrame.
    '''
    tqdm.pandas()
    inputDF.reset_index(inplace=True)

    d = dict(zip(mappingDF.index, mappingDF.iloc[:, 0]))
    inputDF.iloc[:, 0] = inputDF.iloc[:, 0].progress_map(
        lambda x: d.get(x, np.nan))

    inputDF.dropna(inplace=True, subset=[inputDF.columns[0]])
    inputDF.set_index(inputDF.columns[0], inplace=True)


def createTertiaryMatrix(inputDF):
    df = inputDF.copy()

    def mapter(x):
        if x >= 0.95:
            return 1
        elif x <= -0.95:
            return -1
        else:
            return 0

    df = df.applymap(mapter)
    return(df)


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
    if not (lib == 'gene' or lib == 'attribute'):
        return

    filenameGMT = path + name + \
        '_%s.gmt' % str(datetime.date.today())[0:7].replace('-', '_')

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


def createSimilarityMatrix(inputDF, metric):
    '''
    Creates a similarity matrix between the rows of the inputDF based on
    the metric specified. The resulting matrix has both rows and columns labeled
    by the index of inputDF.
    '''
    similarity_matrix = dist.pdist(inputDF, metric)
    similarity_matrix = dist.squareform(similarity_matrix)
    similarity_matrix = 1 - similarity_matrix
    similarity_df = pd.DataFrame(
        data=similarity_matrix, index=inputDF.index, columns=inputDF.index)
    similarity_df.index.name = None
    similarity_df.columns.name = None
    return(similarity_df)


def createGeneList(inputDF):
    '''
    Creates a list of genes and the corresponding Entrez Gene IDs(supplied by
    the NCBI)

    Note: this differs from the previous function in its behavior with dealing
    with genes that do not have an ID. This function will set the id of the gene
    to -1, whereas the previous script will set them to np.nan.
    '''
    gene_list = inputDF.index
    d = dict(zip(getGeneIDs.index, getGeneIDs.iloc[:, 0]))

    gene_ids = np.array([d.get(x, -1) if np.isfinite(d.get(x, -1))
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
    return(attribute_list)


def createGeneAttributeEdgeList(inputDF, attributelist, genelist, path, name):
    '''
    Creates the gene-attribute edge list from the given input DataFrame,
    attribute and gene lists. The year and month are added at the
    end of the name. The path the file is saved to is thus
        path + name + '_<year>_<month>.gmt'
    Also prints the number of cells in inputDF that are statistically
    significant, i.e. > 0.95 confidence.
    Requires:
        attributelist and genelist were generated from running
        createAttributeList and createGeneList on inputDF, respectively.
    '''
    filenameGMT = path + name + \
        '_%s.tsv' % str(datetime.date.today())[0: 7].replace('-', '_')

    temp = pd.DataFrame(columns=['GeneSym', 'GeneID', 'Attribute', 'Weight'])

    with open(filenameGMT, 'w') as f:
        temp['GeneSym'] = [x for x in genelist['GeneSym'].values.tolist()
                           for _ in range(len(attributelist))]
        temp['GeneID'] = [x for x in genelist['GeneID'].values.tolist()
                          for _ in range(len(attributelist))]
        temp['Attribute'] = attributelist.index.values.tolist() * len(genelist)
        temp['Weight'] = inputDF.values.flatten()
        temp.to_csv(f, index=False, sep='\t')

    arr = inputDF.to_numpy()
    count = np.sum(arr >= 0.95) + np.sum(arr <= -0.95)
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
    return(newDF)


def saveData(inputDF, path, name, x=None, y=None, symmetric=False, compression=''):
    pass


def createBinaryMatrix(inputDF, ppi=False):

    if ppi:

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
        genes = list(set(inputDF.iloc[:, 0].unique().tolist()))

        attributes = list(set(inputDF.iloc[:, 1].unique().tolist()))

        matrix = pd.DataFrame(index=genes, columns=attributes, data=0.0)

        for i, gene in enumerate(tqdm(genes)):

            lst = inputDF.loc[(inputDF.iloc[:, 0] == gene),
                              inputDF.columns[1]].values.tolist()

            matrix.at[gene, lst] = 1

        return(matrix)
