# get sparsity and total number of features
# in theory there should be a lot of 0 values in the weights because the sparsity priors but that does not seem to be the case
# we can obtain the sparsity of the original data, measuring the number of zero values relative to the dimensions of the data
# are there NA in our continuous data? no
import os
import pandas as pd

INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
PATHS = {
    "Metabolomics": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv",
    "Lipidomics": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv",
    "Microbiota": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv",
    "SV_DEL": f"{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv",
    "SV_DUP": f"{INPUT_DIR}/SV_DUP_Patient_Gene_Matrix.tsv",
    "SV_INS": f"{INPUT_DIR}/SV_INS_Patient_Gene_Matrix.tsv",
    "SV_INV": f"{INPUT_DIR}/SV_INV_Patient_Gene_Matrix.tsv",
    "SV_TRA": f"{INPUT_DIR}/SV_TRA_Patient_Gene_Matrix.tsv",
    "VC_11": f"{INPUT_DIR}/MPA_GT_1_1.tsv",
    "VC_12": f"{INPUT_DIR}/MPA_GT_1_2.tsv",
    "EXPRESSION": f"{INPUT_DIR}/tpmexpression.tsv"
}
outputpath = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/featuresandsparsity"
os.makedirs(outputpath, exist_ok=True)

numberfeaturesdict = {}
sparsitydict = {}

def obtainnumberfeatures(data):
    numberfeatures = len(data.columns.tolist())
    return numberfeatures

def obtainsparsity(data):
    totalzeros = (data == 0).sum().sum()
    total = data.shape[0]*data.shape[1]
    totalsparsity = (totalzeros / total) * 100
    return totalsparsity

for view, path in PATHS.items():
    data = pd.read_csv(path, sep="\t", index_col=0)

    numberfeatures = obtainnumberfeatures(data)
    numberfeaturesdict[view] = numberfeatures

    totalsparsity = obtainsparsity(data)
    sparsitydict[view] = totalsparsity

numberfeaturesdf = pd.DataFrame.from_dict(numberfeaturesdict, orient="index").reset_index()
numberfeaturesdf.columns = ["View", "Total number of features"]
sparsitydf = pd.DataFrame.from_dict(sparsitydict, orient="index").reset_index()
sparsitydf.columns = ["View", "Sparsity (%)"]
finaldf = pd.merge(numberfeaturesdf, sparsitydf, on="View")
print(finaldf)
finaldf.to_csv(os.path.join(outputpath, "tableuniquevaluessparsity.csv"), index=False)
