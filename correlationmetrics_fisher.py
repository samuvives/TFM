# normalise the correlation matrix using the fisher transformation
# apply the fisher transformation to the correlation data before representing it in an histogram
# get the file from the factors directory
# get the correlation matrix
# apply the fisher transformation to the correlation matrix
# do an histogram of the correlation values
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


path = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS"
OUTPUT_PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/histcorrelationanalysis"
listallfactors = []

def getcorrelations(file):
    data = pd.read_csv(file, index_col = 0)
    correlateddata = data.corr(method="spearman")
    return correlateddata


def fishertransformation(r):
    # z = 1/2 * ln((1+r)/(1-r))
    # revisar si sirve np.arctanh
    transf_r = 0.5 * np.log((1 + r) / (1 - r))
    return transf_r
    

def processcorrmatrix(corr_matrix):
    # replace spaces
    corr_matrix.columns = corr_matrix.columns.str.replace(' ', '_')
    corr_matrix.index = corr_matrix.index.str.replace(' ', '_')
    
    # obtain a table
    mask = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
    datafrommatrix = corr_matrix.where(mask).stack().reset_index()
    datafrommatrix.columns = ['F1', 'F2', 'corrbarvalues']
    datafrommatrix['corrbarnames'] = datafrommatrix['F1'] + " vs " + datafrommatrix['F2']
    return datafrommatrix


def preparingplot(df, namecolumn, numberfactors):
    # mean de todos los valores absolutos en namecolumn
    mean = df[namecolumn].abs().mean()
    median = df[namecolumn].abs().median()
    mode = df[namecolumn].abs().mode()
    numModes = mode.shape[0]
    if numModes > 1:
        print(f"{numberfactors}{namecolumn}{numModes}")
    mode = mode[0]
    listallfactors.append({
        "numberfactors": numberfactors,
        "mean": mean,
        "median": median,
        "mode": mode
        })
    

def corrfactorshist(df, namecolumn, hist_title, save_path):
    """Makes the histograms for each of the runs"""
    plt.figure(figsize=(10, 6))
    sns.histplot(df[namecolumn], binwidth=0.05, kde=True, color='#f26a8d')

    plt.title(hist_title)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def getcorrelationsfromfiles(numberfactors, path=path):
    filepath = os.path.join(path, numberfactors)
    filepath = os.path.join(filepath, f"complete_factors_Z_{numberfactors}.csv")

    cormatrix = getcorrelations(filepath)

    fishercormatrix = fishertransformation(cormatrix)

    datafrommatrix = processcorrmatrix(fishercormatrix)

    OUTPUT_FILE = os.path.join(OUTPUT_PATH, f"corrfisherhist_{numberfactors}.png")
    corrfactorshist(datafrommatrix, "corrbarvalues", f"corrfisherhist_{numberfactors}", OUTPUT_FILE)

    preparingplot(datafrommatrix, "corrbarvalues", numberfactors)

def linescatterplot(column, color):
    fig, ax = plt.subplots()
    ax.scatter(dfallfactors["numberfactors"], dfallfactors[column], color=color)
    ax.plot(dfallfactors["numberfactors"], dfallfactors[column], color=color)
    ax.set_xlabel('Number of factors')

    ylabel = column.title() + " of the correlation" 
    ax.set_ylabel(ylabel)

    title = column.title() + " absolute correlation of the factors through different number of factors\n(Adjusted by Fisher transformation)"
    ax.set_title(title)

    plt.tight_layout()

    save = "fishercorrperfactors.png"
    save = column + save
    save_path = os.path.join(OUTPUT_PATH, save)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


alldirectories = [f for f in os.listdir(path) if f.startswith("K")]
for directory in alldirectories:
    getcorrelationsfromfiles(directory)

# plots across the factors
# order the factors
dfallfactors = pd.DataFrame(listallfactors)
dfallfactors["num"] = dfallfactors["numberfactors"].str.replace("K", "").astype(int)
dfallfactors = dfallfactors.sort_values(by="num").reset_index(drop=True)


# mean
linescatterplot("mean", #ff6b35)

# median
linescatterplot("median", #2a9d8f)

# mode
linescatterplot("mode", #e9c46a)

