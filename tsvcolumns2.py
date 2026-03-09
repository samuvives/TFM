# script to add sufixes in tsv files, applied to metabolomics, lipidomics and microbiota data
import pandas as pd
import os

config = {
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_microbiota_SOFA_case_tumoral.tsv": "MIC",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_Metabolomics_data_case.tsv": "MET",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_lipidomics_data_case.tsv": "LIP"
}

def addsufixes_tsv(filepath, suffix):
    if os.path.exists(filepath):
        try:
            df = pd.read_csv(filepath, sep='\t')
            
            newcolumns = [
                "sample_ID" if i == 0 else f"{col}_{suffix}" 
                for i, col in enumerate(df.columns)
            ]
            
            df.columns = newcolumns
            
            df.to_csv(filepath, sep='\t', index=False)
            print(f"File with path '{filepath}' processed. Added suffix: _{suffix}")
            
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
    else:
        print(f"Error: path '{filepath}' not found")

for filepath, suffix in config.items():
    addsufixes_tsv(filepath, suffix)
