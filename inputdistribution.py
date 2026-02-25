# do not do an histogram, because these are too many values
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
OUTPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/inputdistribution"
os.makedirs(OUTPUT_DIR, exist_ok=True)

VIEWS = {
    "Metabolomics": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv",
    "Lipidomics": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv",
    "Microbiota": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv",
    "EXPRESSION_TPMs": f"{INPUT_DIR}/tpmexpression.tsv"
}


def inputboxplot(view, data, mode="original"):
    fig, ax = plt.subplots()
    # todo
    sns.boxplot(data=data.T, orient="h", showfliers=False)
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
    plt.savefig(os.path.join(OUTPUT_DIR, f"horilollipatients_{view}.png"))


for view, path in VIEWS.items():
    data = pd.read_csv(path, sep="\t", index_col=0)
    inputboxplot(view, data)
    createstatstable(view, data)
    patientsrange(view, data)

statstable = pd.concat(all_stats, axis=1)
statstable.to_csv(os.path.join(OUTPUT_DIR, "inputstatstable.csv"))

# weird patients in transcriptomics
# do a compound histogram
weirdpatientstranscriptomic = ["4772183", "4783769", "5083795", "5270792"]
expressiondata = pd.read_csv(VIEWS["EXPRESSION_TPMs"], sep="\t", index_col=0)
expressiondata.index = expressiondata.index.astype(str)
print(expressiondata.head(10))
print("------")
# pd series per patient
expressiondata1 = expressiondata.loc["4772183"] 
expressiondata2 = expressiondata.loc["4783769"] 
expressiondata3 = expressiondata.loc["5083795"] 
expressiondata4 = expressiondata.loc["5270792"] 

fig, axes = plt.subplots(2, 2)

sns.histplot(expressiondata1, kde=True, color="#ffca3a", ax=axes[0, 0])
axes[0, 0].set_title("4772183")

sns.histplot(expressiondata2, kde=True, color="#8ac926", ax=axes[0, 1])
axes[0, 1].set_title("4783769")

sns.histplot(expressiondata3, kde=True, color="#1982c4", ax=axes[1, 0])
axes[1, 0].set_title("5083795")

sns.histplot(expressiondata4, kde=True, color="#6a4c93", ax=axes[1, 1])
axes[1, 1].set_title("5270792")

for ax in axes.flatten():
    ax.set_xlabel("Values")
    ax.set_ylabel("Frequency")

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "weirdexpressionpatients.png"))
