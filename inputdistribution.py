# do not do an histogram, because these are too many values
import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
OUTPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/inputdistribution"
os.makedirs(OUTPUT_DIR, exist_ok=True)

CONTINUOUSVIEWS = {
    "Metabolomics": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv",
    "Lipidomics": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv",
    "Microbiota": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv",
    "EXPRESSION_TPMs": f"{INPUT_DIR}/tpmexpression.tsv"
}

def inputboxplot(view, data, mode="original"):
    fig, ax = plt.subplots()
    if mode == "withoutoutliers":
        sns.boxplot(data=data.T, orient="h", showfliers=False, fliersize=2)
    else:
        sns.boxplot(data=data.T, orient="h", showfliers=True, fliersize=2)
    ax.set_title(f"Data of {mode} {view} per patient")
    ax.set_xlabel("Values")
    ax.set_ylabel("Patient")
    ax.tick_params(axis="y", labelsize=5)
    plt.tight_layout()
    OUTPUTPATH = os.path.join(OUTPUT_DIR, f"boxplot_{mode}_{view}.png")
    plt.savefig(OUTPUTPATH)
    plt.close()


# obtain the 
# to do the horizontal lolli with the min and max
# to create a table of stats of the input

all_stats = []

def createstatstable(view, data):
    stats = pd.DataFrame({
        f"{view}_Min": data.min(axis=1),
        f"{view}_Max": data.max(axis=1),
        f"{view}_Mean": data.mean(axis=1),
        f"{view}_Median": data.median(axis=1),
        f"{view}_Std": data.std(axis=1),
        f"{view}_Zeros": (data == 0).sum(axis=1),
        f"{view}_ZerosProportion": ((data == 0).sum(axis=1)) / data.shape[1],
        f"{view}_NAs": (data.isnull().sum(axis=1)) / data.shape[1],
        f"{view}_NAsProportion": (data.isnull().sum(axis=1)) / data.shape[1]
    })
    all_stats.append(stats)


def patientsrange(view, data):
    """Poner minimo maximo"""
    datastats = pd.concat([data.min(axis=1), data.max(axis=1)], axis=1)
    datastats.columns = ["Min", "Max"]
    numpatients = range(datastats.shape[0])
    fig, ax = plt.subplots()
    ax.scatter(datastats["Min"], numpatients, color="Black")
    ax.scatter(datastats["Max"], numpatients, color="Black")
    ax.hlines(y=numpatients, xmin=datastats["Min"], xmax=datastats["Max"], color="Black")
    ax.set_yticks(numpatients)
    ax.set_yticklabels(datastats.index)
    ax.set_ylabel("Patients")
    ax.tick_params(axis="y", labelsize=5)
    ax.set_xlabel("Min-Max value")
    ax.set_title(f"Min and Max values of the view {view}")
    plt.savefig(os.path.join(OUTPUT_DIR, f"horilollipatients_{view}.png"))


for view, path in CONTINUOUSVIEWS.items():
    data = pd.read_csv(path, sep="\t", index_col=0)
    inputboxplot(view, data)
    inputboxplot(view, data, mode="withoutoutliers")

    createstatstable(view, data)
    patientsrange(view, data)
    scaleddata = pd.DataFrame(StandardScaler().fit_transform(data), index=data.index, columns=data.columns)
    inputboxplot(view, scaleddata, mode="scaled")


statstable = pd.concat(all_stats, axis=1)
statstable.to_csv(os.path.join(OUTPUT_DIR, "inputstatstable.csv"))

# weird patients in transcriptomics
# do a compound histogram
weirdpatientstranscriptomic = ["4772183", "4783769", "5083795", "5270792"]
expressiondata = pd.read_csv(CONTINUOUSVIEWS["EXPRESSION_TPMs"], sep="\t", index_col=0)
expressiondata.index = expressiondata.index.astype(str)
print(expressiondata.head(10))
print("------")
# pd series per patient
expressiondata1 = expressiondata.loc["4772183"] 
expressiondata2 = expressiondata.loc["4783769"] 
expressiondata3 = expressiondata.loc["5083795"] 
expressiondata4 = expressiondata.loc["5270792"] 

mapping_df = pd.read_excel("/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/CORRESPONDENCIA.xlsx")
mapping_dict = dict(zip(mapping_df["NEWEXPRESSION"], mapping_df["EXPRESSION"]))
mapping_dict = {str(k): v for k, v in mapping_dict.items()}

# all again
expressiondatamapped = expressiondata.rename(index=mapping_dict)
inputboxplot("allpatientsmapped", expressiondatamapped)
inputboxplot("allpatientsmapped", expressiondatamapped, mode="withoutoutliers")
fig, ax = plt.subplots()
sns.boxplot(data=expressiondatamapped.T, orient="h", showfliers=False, fliersize=2)
ax.set_title(f"EXPRESSION data per patient")
ax.set_xlabel("Values")
ax.set_ylabel("Patient")
ax.tick_params(axis="y", labelsize=5)
plt.tight_layout()
OUTPUTTOSEND = os.path.join(OUTPUT_DIR, f"boxplot_EXPRESSION_tosend.png")
plt.savefig(OUTPUTTOSEND)
plt.close()

# only weird patients
weirdpatientsboxplotdf = expressiondata.loc[["4772183", "4783769", "5083795", "5270792"]]
weirdpatientsboxplotdf = weirdpatientsboxplotdf.rename(index=mapping_dict)
print(weirdpatientsboxplotdf.index)
inputboxplot("weirdpatientsmapped", weirdpatientsboxplotdf)
inputboxplot("weirdpatientsmapped", weirdpatientsboxplotdf, mode="withoutoutliers")

# filter outliers for histogram
Q1_1 = expressiondata1.quantile(0.25)
Q3_1 = expressiondata1.quantile(0.75)
IQR_1 = Q3_1 - Q1_1
limite_inferior_1 = Q1_1 - 1.5 * IQR_1
limite_superior_1 = Q3_1 + 1.5 * IQR_1
expressiondata1 = expressiondata1[(expressiondata1 >= limite_inferior_1) & (expressiondata1 <= limite_superior_1)]

Q1_2 = expressiondata2.quantile(0.25)
Q3_2 = expressiondata2.quantile(0.75)
IQR_2 = Q3_2 - Q1_2
limite_inferior_2 = Q1_2 - 1.5 * IQR_2
limite_superior_2 = Q3_2 + 1.5 * IQR_2
expressiondata2 = expressiondata2[(expressiondata2 >= limite_inferior_2) & (expressiondata2 <= limite_superior_2)]

Q1_3 = expressiondata3.quantile(0.25)
Q3_3 = expressiondata3.quantile(0.75)
IQR_3 = Q3_3 - Q1_3
limite_inferior_3 = Q1_3 - 1.5 * IQR_3
limite_superior_3 = Q3_3 + 1.5 * IQR_3
expressiondata3 = expressiondata3[(expressiondata3 >= limite_inferior_3) & (expressiondata3 <= limite_superior_3)]

Q1_4 = expressiondata4.quantile(0.25)
Q3_4 = expressiondata4.quantile(0.75)
IQR_4 = Q3_4 - Q1_4
limite_inferior_4 = Q1_4 - 1.5 * IQR_4
limite_superior_4 = Q3_4 + 1.5 * IQR_4
expressiondata4 = expressiondata4[(expressiondata4 >= limite_inferior_4) & (expressiondata4 <= limite_superior_4)]

fig, axes = plt.subplots(2, 2)

sns.histplot(expressiondata1, kde=True, color="#ffca3a", ax=axes[0, 0], bins=100)
axes[0, 0].set_yscale('log')
axes[0, 0].set_title("4772183")

sns.histplot(expressiondata2, kde=True, color="#8ac926", ax=axes[0, 1], bins=100)
axes[0, 1].set_yscale('log')
axes[0, 1].set_title("4783769")

sns.histplot(expressiondata3, kde=True, color="#1982c4", ax=axes[1, 0], bins=100)
axes[1, 0].set_yscale('log')
axes[1, 0].set_title("5083795")

sns.histplot(expressiondata4, kde=True, color="#6a4c93", ax=axes[1, 1], bins=100)
axes[1, 1].set_yscale('log')
axes[1, 1].set_title("5270792")

for ax in axes.flatten():
    ax.set_xlabel("Values")
    ax.set_ylabel("Frequency")

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "weirdexpressionpatients.png"))

#-----------------------------------------------------------------------
# BINARY DATA
#-----------------------------------------------------------------------

def scatterbinary(data, view):
    mutatedproportion = data.sum(axis=0) / data.shape[0]
    fig, ax = plt.subplots()
    ax.scatter(x=range(len(mutatedproportion)), 
        y=mutatedproportion.values, s=0.5)
    ax.set_xlabel("Features")
    ax.set_xticks([])
    ax.set_xticklabels([])

    ax.set_ylabel("Relative frequency of mutation in patients")
    ax.set_title(f"Frequency of mutation in the patients. {view}")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"scatterbinary_{view}.png"))

BINARYVIEWS = {
    "SV_DEL": f"{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv",
    "SV_DUP": f"{INPUT_DIR}/SV_DUP_Patient_Gene_Matrix.tsv",
    "SV_INS": f"{INPUT_DIR}/SV_INS_Patient_Gene_Matrix.tsv",
    "SV_INV": f"{INPUT_DIR}/SV_INV_Patient_Gene_Matrix.tsv",
    "SV_TRA": f"{INPUT_DIR}/SV_TRA_Patient_Gene_Matrix.tsv",
    "VC_11": f"{INPUT_DIR}/MPA_GT_1_1.tsv",
    "VC_12": f"{INPUT_DIR}/MPA_GT_1_2.tsv",
}


for view, path in BINARYVIEWS.items():
    data = pd.read_csv(path, sep="\t", index_col=0)
    scatterbinary(data, view)

# count the 1s per column / numberpatients
