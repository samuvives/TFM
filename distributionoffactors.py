# check how the values of the Z matrix are distributed
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

path = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K12"
inputpath = os.path.join(path, "complete_factors_Z_K12.csv")
outputdir = os.path.join(path, "factorsdistribution")
os.makedirs(outputdir, exist_ok=True)

# load the Z matrix
df = pd.read_csv(inputpath, index_col=0)

# function to make an histogram
def factorhist(factor, colname):
    colname = colname.replace(" ", "_")
    fig, ax = plt.subplots()
    sns.histplot(
        data=factor,
        binwidth=0.05,
        kde=True,
        ax=ax,
        color="#7f5539")
    ax.set_title(f"Histogram of the Z matrix. K12. {colname}")
    ax.set_xlabel("Values")
    ax.set_ylabel("Frequency")
    plt.tight_layout()
    outputpath = os.path.join(outputdir, f"disthist_{colname}.png")
    plt.savefig(outputpath)
    plt.close()
    

# apply the function to each column of Z (factor)
for colname in df.columns:
    factor = df[colname]
    factorhist(factor, colname)
