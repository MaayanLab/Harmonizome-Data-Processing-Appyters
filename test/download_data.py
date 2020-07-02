# Script to acquire RNA-Seq and clinical data from TCGA
import requests
import json
import os
import re
import gzip
import shutil
import tarfile
import pathlib
import pandas as pd
import numpy as np
from maayanlab_bioinformatics.harmonization import ncbi_genes
import math

# Endpoints
base_url = 'https://api.gdc.cancer.gov/'
files_endpt = base_url + 'files/'
genes_endpt = base_url + 'genes/'
cases_endpt = base_url + 'cases/'
data_endpt = base_url + "data/"

# data type of files we want
data_type = "htseq.counts"

# minimum number of cases for a single cancer
min_cases = 150

# cancer types listed under the primary_diagnosis field in TCGA with at least one RNA seq file entry
cancer_types = [
  'Pleomorphic carcinoma',
  'Adrenal cortical carcinoma',
  'Amelanotic melanoma',
  'Basaloid squamous cell carcinoma',
  'Adenosquamous carcinoma',
  'Squamous cell carcinoma, small cell, nonkeratinizing',
  'Papillary adenocarcinoma, NOS',
  'Superficial spreading melanoma',
  'Pheochromocytoma, malignant',
  'Mucinous adenocarcinoma, endocervical type',
  'Thymoma, type B1, NOS',
  'Neuroendocrine carcinoma, NOS',
  'Epithelioid mesothelioma, malignant',
  'Thymoma, type B2, malignant',
  'Medullary carcinoma, NOS',
  'Paget disease and infiltrating duct carcinoma of breast',
  'Oligodendroglioma, anaplastic',
  'Squamous cell carcinoma, keratinizing, NOS',
  'Phyllodes tumor, malignant',
  'Intraductal micropapillary carcinoma',
  'Leiomyosarcoma, NOS',
  'Apocrine adenocarcinoma',
  'Large cell neuroendocrine carcinoma',
  'Papillary serous cystadenocarcinoma',
  'Thymoma, type B3, malignant',
  'Solid carcinoma, NOS',
  'Oligodendroglioma, NOS',
  'Mixed epithelioid and spindle cell melanoma',
  'Thymoma, type A, NOS',
  'Adenocarcinoma with mixed subtypes',
  'Renal cell carcinoma, chromophobe type',
  'Infiltrating duct mixed with other types of carcinoma',
  'Aggressive fibromatosis',
  'Papillary carcinoma, NOS',
  'Bronchiolo-alveolar carcinoma, non-mucinous',
  'Infiltrating duct and lobular carcinoma',
  'Pleomorphic liposarcoma',
  'Carcinoma, NOS',
  'Giant cell sarcoma',
  'Cribriform carcinoma, NOS',
  'Yolk sac tumor',
  'Abdominal fibromatosis',
  'Thymoma, type AB, NOS',
  'Dedifferentiated liposarcoma',
  'Teratoma, malignant, NOS',
  'Signet ring cell carcinoma',
  'Pheochromocytoma, NOS',
  'Synovial sarcoma, biphasic',
  'Fibrous mesothelioma, malignant',
  'Papillary squamous cell carcinoma',
  'Infiltrating duct carcinoma, NOS',
  'Myxoid leiomyosarcoma',
  'Epithelioid cell melanoma',
  'Endometrioid adenocarcinoma, NOS',
  'Cholangiocarcinoma',
  'Extra-adrenal paraganglioma, malignant',
  'Nonencapsulated sclerosing carcinoma',
  'Thymoma, type B1, malignant',
  'Papillary carcinoma, columnar cell',
  'Extra-adrenal paraganglioma, NOS',
  'Paraganglioma, NOS',
  'Mixed glioma',
  'Mesothelioma, malignant',
  'Spindle cell melanoma, NOS',
  'Malignant melanoma, NOS',
  'Seminoma, NOS',
  'Liposarcoma, well differentiated',
  'Follicular adenocarcinoma, NOS',
  'Follicular carcinoma, minimally invasive',
  'Squamous cell carcinoma, large cell, nonkeratinizing, NOS',
  'Papillary transitional cell carcinoma',
  'Spindle cell melanoma, type B',
  'Adenocarcinoma, intestinal type',
  'Squamous cell carcinoma, NOS',
  'Acral lentiginous melanoma, malignant',
  'Transitional cell carcinoma',
  'Papillary carcinoma, follicular variant',
  'Bronchio-alveolar carcinoma, mucinous',
  'Nodular melanoma',
  'Renal cell carcinoma, NOS',
  'Teratocarcinoma',
  'Astrocytoma, anaplastic',
  'Undifferentiated sarcoma',
  'Intraductal papillary adenocarcinoma with invasion',
  'Bronchiolo-alveolar adenocarcinoma, NOS',
  'Acute myeloid leukemia, NOS',
  'Thymoma, type A, malignant',
  'Tubular adenocarcinoma',
  'Serous cystadenocarcinoma, NOS',
  'Desmoplastic melanoma, malignant',
  'Carcinosarcoma, NOS',
  'Thymoma, type B2, NOS',
  'Oxyphilic adenocarcinoma',
  'Adenocarcinoma, endocervical type',
  'Hepatocellular carcinoma, fibrolamellar',
  'Adenocarcinoma with neuroendocrine differentiation',
  'Mesodermal mixed tumor',
  'Mesothelioma, biphasic, malignant',
  'Clear cell adenocarcinoma, NOS',
  'Astrocytoma, NOS',
  'Glioblastoma',
  'Adenoid cystic carcinoma',
  'Micropapillary carcinoma, NOS',
  'Synovial sarcoma, NOS',
  'Adenocarcinoma, NOS',
  'Infiltrating lobular mixed with other types of carcinoma',
  'Embryonal carcinoma, NOS',
  'Carcinoma, diffuse type',
  'Hepatocellular carcinoma, clear cell type',
  'Squamous cell carcinoma, spindle cell',
  'Synovial sarcoma, spindle cell',
  'Combined hepatocellular carcinoma and cholangiocarcinoma',
  'Endometrioid adenocarcinoma, secretory variant',
  'Malignant peripheral nerve sheath tumor',
  'Thymoma, type AB, malignant',
  'Metaplastic carcinoma, NOS',
  'Fibromyxosarcoma',
  'Thymic carcinoma, NOS',
  'Mullerian mixed tumor',
  'Teratoma, benign',
  'Carcinoma, undifferentiated, NOS',
  'Adenocarcinoma in tubolovillous adenoma',
  'Paraganglioma, malignant',
  'Acinar cell carcinoma',
  'Mucinous adenocarcinoma',
  'Secretory carcinoma of breast',
  'Lentigo maligna melanoma',
  'Malignant lymphoma, large B-cell, diffuse, NOS',
  'Serous surface papillary carcinoma',
  'Basal cell carcinoma, NOS',
  'Malignant fibrous histiocytoma',
  'Hepatocellular carcinoma, NOS',
  'Lobular carcinoma, NOS',
  'Mixed germ cell tumor'
]

os.makedirs("TCGA_cancer_downloads", exist_ok = True)

'''
Loop over all cancer types, creating the following files:
- /data
    - /{cancer_type}
        - all RNA-seq files for this cancer type
        - data.csv : single file with all RNA seq data for this cancer type
        - clinical_data.csv : single file with all clinical data for this cancer type
- /TCGA_cancer_downloads
    - all .tar.gz files
    - all unzipped file folders
        - each containing a single RNA-seq file
'''

for cancer_type in cancer_types:
    ############# RNA-seq DATA ################

    # Retrieve data from the TCGA API

    # data type of files we want
    data_type = "htseq.counts"

    # The 'fields' parameter is passed as a comma-separated string of single names.
    fields = "file_id,file_name,cases.case_id"

    # filter files for only RNA-Seq results
    filters = {
        "op": "and",
         "content":[
             {
                "op": "in",
                "content":
                 {
                     "field": "files.experimental_strategy",
                     "value": ["RNA-Seq"],
                 }
             },
             {
                "op": "in",
                "content":
                 {
                     "field": "access",
                     "value": ["open"],

                 }
             },
             {
                "op": "in",
                "content":
                 {
                     "field": "files.file_name",
                     "value": ["*htseq.counts.gz"],
                 }
             },

             {
                "op": "in",
                "content":
                 {
                     "field": "cases.diagnoses.primary_diagnosis",
                     "value": [cancer_type],
                 }
             },
         ],
    }

    # build parameters object
    params = {
        "fields": fields,
        "filters": json.dumps(filters),
        "size": 100000
    }

    # get list of all files with RNA-seq results
    response = requests.get(files_endpt, params = params) # optionally also provide params argument
    data = json.loads(response.content.decode("utf-8"))

    # get list of results
    results = data["data"]["hits"]
    results = list(filter(lambda x: data_type in x["file_name"], results))

    file_uuid_list = [ entry["file_id"] for entry in results]
    case_uuid_list = [ entry["cases"][0]["case_id"] for entry in results]

    num_cases = len(case_uuid_list)
    print(f"{cancer_type}: {num_cases} cases")

    # continue to next cancer type if this case list is too small
    if (num_cases < min_cases): continue

    os.makedirs(f'data/{cancer_type}', exist_ok = True)


    df_files_cases=pd.DataFrame({"case": case_uuid_list },  index=file_uuid_list)
    file_to_case = df_files_cases.to_dict()["case"] # dict mapping file id to case id

    params = {"ids": file_uuid_list}

    # A POST is used, so the filter parameters can be passed directly as a Dict object.
    response = requests.post(data_endpt,
                            data = json.dumps(params),
                            headers={
                                "Content-Type": "application/json"})

    # filename is found in the Content-Disposition header of response
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    downloads_folder = "TCGA_cancer_downloads/"

    # Save .tar.gz zipped file to TCGA_downloads folder
    with open(downloads_folder + file_name, "wb") as f_out:
        f_out.write(response.content)

    # extract the root tar archive
    tar = tarfile.open(downloads_folder + file_name, "r:gz")
    tar.extractall(f'./{downloads_folder}')
    folder = file_name.split(".tar.gz")[0]

    for tarinfo in tar:
        if (tarinfo.name == "MANIFEST.txt"): continue
        file_id = tarinfo.name.split("/")[0]

        # unzip inner .gz files
        with gzip.open(downloads_folder + tarinfo.name, "rb") as f_in:
            os.makedirs(f"data/{cancer_type}",exist_ok = True)
            with open(f"data/{cancer_type}/{file_to_case[file_id]}.txt", "wb") as f_out:
                f_out.write(f_in.read())

    tar.close()

    # initialize empty df
    df = pd.DataFrame({"gene": []})
    df = df.set_index("gene")

    # loop over files, merging with pre-existing data
    for file in pathlib.Path(f'data/{cancer_type}').glob('*.txt'):
        with open(file, "rb") as f_in:
            new_df = pd.read_csv(f_in, sep = "\t", header = None)
            file_id = re.findall(f"data/{cancer_type}/(.+).txt", f_in.name)[0]
            new_df.columns = ["gene", file_id]
            new_df.gene.replace(to_replace = r'\..*$', value = "", regex=True,
               inplace=True) # collapse all versions of same gene to one gene
            new_df = new_df.set_index("gene")
            df = pd.DataFrame.merge(df, new_df, how="outer", left_on = "gene", right_on = "gene")

    # drop rows not corresponding to genes (i.e. metadata)
    non_genes = list(filter(lambda val: not "ENSG" in val,np.array(df.index.values)))
    df = df.drop(non_genes)

    ncbi = pd.DataFrame(ncbi_genes.ncbi_genes_fetch())

    def get_ensembl_id(ids):
        ids = "".join(ids)
        ensembl = re.findall("Ensembl:(.*)", ids)
        if (len(ensembl) == 1):
            return ensembl[0]
        else:
            return None

    all_ids = ncbi.dbXrefs.values
    ensembl_ids = [ get_ensembl_id(ids) for ids in all_ids]

    ncbi = ncbi[["dbXrefs", "Symbol", "type_of_gene"]]
    ncbi["ensembl"] = ensembl_ids
    ncbi = ncbi.drop(columns=["dbXrefs"])
    ncbi = ncbi.set_index("ensembl")

    ensembl_to_symbol = ncbi.to_dict()["Symbol"]
    ensembl_to_gene_type = ncbi.to_dict()["type_of_gene"]


    data_ensembl_ids = df.index.to_list()

    # if the key is present, return it; otherwise, set the index for the corresponding row as its ensembl id
    def id_to_symbol(key):
        if (key in ensembl_to_symbol):
            return ensembl_to_symbol[key]
        else:
            return key # can return key here to maintain some gene identity

    def id_to_type(key):
        if (key in ensembl_to_gene_type):
            return ensembl_to_gene_type[key]
        else:
            return None

    data_symbols = [ id_to_symbol(key) for key in data_ensembl_ids ]
    data_types = [ id_to_type(key) for key in data_ensembl_ids ]

    df["symbol"] = data_symbols
    df["type_of_gene"] = data_types
    df_symbol = df.set_index("symbol")

    # drop non protein-coding genes
    df_symbol = df_symbol[df_symbol["type_of_gene"] == "protein-coding"]
    df_symbol = df_symbol.drop(columns=["type_of_gene"])
    df_symbol.head()

    # write the final csv, ready for normalization and further preprocessing
    df_symbol.to_csv(f'data/{cancer_type}/data.csv', encoding='utf-8')


    ############# Clinical Data ################

    # get all demographic and diagnoses fields
    cases_fields = requests.get(cases_endpt + "_mapping").json()["fields"]
    demographic_fields = list(filter(lambda x: "demographic" in x, cases_fields))
    diagnoses_fields = list(filter(lambda x: "diagnoses" in x, cases_fields))


    # initialize dataframe
    demographic_column_names = [ field.split(".")[1] for field in demographic_fields ]
    diagnoses_column_names = [ field.split(".")[1] for field in diagnoses_fields ]

    columns = list(set([*demographic_column_names,*diagnoses_column_names]))
    df_clinical = pd.DataFrame({}, columns=columns)
    df_clinical["case_id"] = []

    # get demographics and diagnosis data for each case,
    # merging with pre-exisiting dataframe
    for case in case_uuid_list:
        fields=",".join([*demographic_fields, *diagnoses_fields])
        params={
            "fields": fields
        }
        response = requests.get(cases_endpt + case, params=params).json()["data"]
        demographic_data = response["demographic"]
        diagnoses_data = response["diagnoses"]
        diagnoses_data = diagnoses_data[0]
        if "treatments" in diagnoses_data:
            del diagnoses_data["treatments"] # do not load treatment data
        df_case = pd.DataFrame({**demographic_data,**diagnoses_data}, index=[case])
        df_case.head()
        df_case["case_id"] = case
        df_clinical = pd.concat([df_clinical, df_case], join="outer")

    df_clinical = df_clinical.set_index("case_id")

    # make first column "primary_diagnosis"
    cols = ['primary_diagnosis']  + [col for col in df_clinical.columns.values if col != 'primary_diagnosis']
    df_clinical = df_clinical[cols]

    # save final .csv
    df_clinical.to_csv(f"data/{cancer_type}/clinical_data.csv", encoding='utf-8')
