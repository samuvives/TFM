import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/OBTAININGELBO/MOFAFLEX_FINAL_ANALYSIS/K12/complete_weights"
outputpath = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/OBTAININGELBO/postanalysis/weightanalysis"
os.makedirs(outputpath, exist_ok=True)

WEIGHTSDICT = {
    "Lipidomics": "complete_weights_Lipidomics_K12.csv",
    "SV_INS": "complete_weights_SV_INS_K12.csv",  
    "SV_DEL": "complete_weights_SV_DEL_K12.csv",
    "SV_INV": "complete_weights_SV_INV_K12.csv", 
    "SV_DUP": "complete_weights_SV_DUP_K12.csv",
    "SV_TRA": "complete_weights_SV_TRA_K12.csv",
    "VC_11": "complete_weights_VC_11_K12.csv",
    "VC_12": "complete_weights_VC_12_K12.csv",
    "EXPRESSION": "complete_weights_EXPRESSION_K12.csv",
    "Metabolomics": "complete_weights_Metabolomics_K12.csv",
    "Microbiota": "complete_weights_Microbiota_K12.csv" 
}

for view, weightfile in WEIGHTSDICT.items():
    weighfile = os.path.join(path, weightfile)
    weightdata = pd.read_csv(weightfile, index_col=0)
    weightdata = weightdata.melt()

    fig, ax = plt.subplots()
    sns.histplot(
        data=weightdata,
        kde=True,
        ax=ax,
        color="#5f0f40")
    ax.set_title(f"Histogram of the W matrix. K12. {view}")
    ax.set_xlabel("Values")
    ax.set_ylabel("Frequency")
    plt.tight_layout()
    savepath = os.path.join(outputpath, f"weightsdist_{view}.png")
    plt.savefig(savepath)
    plt.close()
