# compare the weights names that have remained
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

approach = "OBTAININGELBO"
PROJECTDIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/"
INPUTDIR = os.path.join(PROJECTDIR, "MOFAINPUT")
APPROACHDIR = os.path.join(PROJECTDIR, approach)
ANALYSISDIR = os.path.join(APPROACHDIR, "MOFAFLEX_FINAL_ANALYSIS")
FACTORSDIR = [f for f in os.listdir(ANALYSISDIR) if f.startswith("K") and os.path.isdir(os.path.join(ANALYSISDIR, f))]
OUTPUTFILE = os.path.join(APPROACHDIR, "postanalysis/featuresproportion.png")

# format dictinputweights:
# MPA_GT_1_1.tsv: complete_weights_VC_11_K12.csv, 

# for each run compare input features with weight matrix features
# list of rows
alldata = []
for k in FACTORSDIR:
    dictinputweights = {
        "MPA_GT_1_1.tsv": f"complete_weights_VC_11_{k}.csv", 
        "MPA_GT_1_2.tsv": f"complete_weights_VC_12_{k}.csv", 
        "SV_DEL_Patient_Gene_Matrix.tsv": f"complete_weights_SV_DEL_{k}.csv", 
        "SV_DUP_Patient_Gene_Matrix.tsv": f"complete_weights_SV_DUP_{k}.csv", 
        "SV_INS_Patient_Gene_Matrix.tsv": f"complete_weights_SV_INS_{k}.csv", 
        "SV_INV_Patient_Gene_Matrix.tsv": f"complete_weights_SV_INV_{k}.csv", 
        "SV_TRA_Patient_Gene_Matrix.tsv": f"complete_weights_SV_TRA_{k}.csv", 
        "tpmexpression.tsv": f"complete_weights_EXPRESSION_{k}.csv", 
        "renamed_microbiota_SOFA_case_tumoral.tsv": f"complete_weights_Microbiota_{k}.csv", 
        "renamed_Metabolomics_data_case.tsv": f"complete_weights_Metabolomics_{k}.csv", 
        "renamed_lipidomics_data_case.tsv": f"complete_weights_Lipidomics_{k}.csv"
    }
    WEIGHTDIR =  os.path.join(ANALYSISDIR, f"{k}/complete_weights")

    for input, weight in dictinputweights.items():
        pathinput = os.path.join(INPUTDIR, input)
        inputdf = pd.read_csv(pathinput, sep="\t", index_col=0)

        pathweight = os.path.join(WEIGHTDIR, weight)
        weightdf = pd.read_csv(pathweight, index_col=0)

        # obtain the nameview extracting from the name of the weight file
        nameview = weight.replace("complete_weights_", "").replace(".csv", "")
        nameview = re.sub(r'_K\d+$', '', nameview)

        # proportion of features that have remained in that view
        remaining_feat_prop = (inputdf.shape[1] - weightdf.shape[0]) / inputdf.shape[1]

        # create a row for the data and append it to the list
        alldata.append({
                "Run": k,
                "View": nameview,
                "Proportion": remaining_feat_prop
            })

# create a df of the list of rows
proportiondf = pd.DataFrame(alldata)
# column to order the factors
proportiondf["RunNumber"] = proportiondf["Run"].str.replace("K", "").astype(int)

# lineplot for each view across the runs
proportiondf = proportiondf.sort_values(["View", "RunNumber"])
fig, ax = plt.subplots()
sns.lineplot(data=proportiondf, x="RunNumber", y="Proportion", hue='View', 
            marker='s', linewidth=2, ax=ax)
ax.set_title("Proportion of features mantained \n relative to input across runs", fontsize=14)
ax.set_ylabel('Proportion')
ax.set_xlabel('Run')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='View')
plt.tight_layout()
plt.savefig(OUTPUTFILE)
plt.close()
