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
