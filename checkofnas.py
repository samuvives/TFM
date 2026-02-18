# script to check for NAs in different files
import os
import pandas as pd

def checknas(PATH, file):
    """Reads the file and prints the total NAs"""
    fullfile = os.path.join(PATH, file)
    data = pd.read_csv(fullfile, sep="\t")
    numberna = int(data.isnull().sum().sum())
    if numberna > 0:
        print(f"{file} has nas")
    else:
        print(f"{file} has no nas")
    return numberna


# files of raw transcript data
PATH = "/home/vant/Escritorio/TFM/datostfm/tpmssplicing/"
listfiles = [
    "FAMCOLON_16.rsem.merged.transcript_tpm.tsv", "FAMCOLON_19.rsem.merged.transcript_tpm.tsv", 
    "FAMCOLON_21.rsem.merged.transcript_tpm.tsv", "FAMCOLON_27.rsem.merged.transcript_tpm.tsv", 
    "transcript_tpm_matrix_suppa.tsv"
    ]

for file in listfiles:
    numberna = checknas(PATH, file)
    print(f"{file} has {numberna} NAs")


# input files of continuous views of MOFA-FLEX
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
VIEWS = {
    "Metabolomics": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv",
    "Lipidomics": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv",
    "Microbiota": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv",
    "EXPRESSION": f"{INPUT_DIR}/tpmexpression.tsv"
}

for view, file in VIEWS.items():
    numberna = checknas(INPUT_DIR, file)
    print(f"{view} file has {numberna} NAs")
