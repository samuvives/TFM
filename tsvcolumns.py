# Adds the sufixes to each omic file

import pandas as pd
import os

config = {
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_microbiota_SOFA_case_tumoral.tsv": "MIC",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_Metabolomics_data_case.tsv": "MET",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_lipidomics_data_case.tsv": "LIP"
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DEL_Patient_Gene_Matrix.tsv": "DEL",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DUP_Patient_Gene_Matrix.tsv": "DUP",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INS_Patient_Gene_Matrix.tsv": "INS",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INV_Patient_Gene_Matrix.tsv": "INV",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_TRA_Patient_Gene_Matrix.tsv": "TRA"
}

for file, suffix in config.items():
    if os.path.exists(archivo):
        try:
            df = pd.read_csv(file, sep='\t')
            
            # ignore the first column sample_ID
            new_cols = [
                "sample_ID" if i == 0 else f"{col}_{suffix}" for i, col in enumerate(df.columns)
            ]
            
            df.columns = new_cols
            
            df.to_csv(archivo, sep='\t', index=False)
            print(f"File: {file} (Suffix: _{suffix}) processed")
            
        except Exception as e:
            print(f"Error processing {file}: {e}")
    else:
        print(f"The file '{file}' was not found in the directory.")
