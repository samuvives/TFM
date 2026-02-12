# script to check for NAs in different files
import os
import pandas as pd
PATH = "/home/vant/Escritorio/TFM/datostfm/tpmssplicing/"
listfiles = [
    "FAMCOLON_16.rsem.merged.transcript_tpm.tsv", "FAMCOLON_19.rsem.merged.transcript_tpm.tsv", 
    "FAMCOLON_21.rsem.merged.transcript_tpm.tsv", "FAMCOLON_27.rsem.merged.transcript_tpm.tsv", 
    "transcript_tpm_matrix_suppa.tsv"
    ]

def checknas(file):
    fullfile = os.path.join(PATH, file)
    data = pd.read_csv(fullfile, sep="\t")
    numberna = int(data.isnull().sum().sum())
    print(numberna)
    if numberna > 0:
        print(f"{file} has nas")
    else:
        print(f"{file} has no nas")
        

for file in listfiles:
    checknas(file)



import os
import pandas as pd
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
VIEWS = {
    "Metabolomics": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv",
    "Lipidomics": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv",
    "Microbiota": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv",
    "EXPRESSION": f"{INPUT_DIR}/tpmexpression.tsv"
}


nadict = {}
for view, path in VIEWS.items():
    data = pd.read_csv(path, sep="\t")
    numbernas = data.isnull().sum().sum()
    nadict[view] = numbernas

print("Results:")
for view, nanumber in nadict.items():
    print(f"View: {view}, Total NAs: {nanumber}")
