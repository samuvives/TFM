# normalice the correlation matrix using the fisher transformation
# apply the fisher transformation to the correlation data before representing it in an histogram
# get the file from the factors directory
# get the correlation matrix
# apply the fisher transformation to the correlation matrix
# do an histogram of the correlation values
# ahora hazlo modular
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
    corr_matrix.columns = corr_matrix.columns.str.replace(' ', '_')
    corr_matrix.index = corr_matrix.index.str.replace(' ', '_')

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


alldirectories = [f for f in os.listdir(path) if f.startswith("K")]
for directory in alldirectories:
    getcorrelationsfromfiles(directory)

# plots across the factors
# order the factors
dfallfactors = pd.DataFrame(listallfactors)
dfallfactors["num"] = dfallfactors["numberfactors"].str.replace("K", "").astype(int)
dfallfactors = dfallfactors.sort_values(by="num")

# function of the graph
def linescatterplot(column, color, title, ylabel, save):
    fig, ax = plt.subplots()
    ax.scatter(dfallfactors["numberfactors"], dfallfactors[column], color=color)
    ax.plot(dfallfactors["numberfactors"], dfallfactors[column], color=color)
    ax.set_title(title)
    ax.set_xlabel('Number of factors')
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    save_path = os.path.join(OUTPUT_PATH, save)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# mean
meantitle = "Mean absolute correlation of the factors through different number of factors\n(Adjusted by Fisher transformation)"
meansave = "meanfishercorrperfactors.png"
linescatterplot("mean", #ff6b35, meantitle,"Mean correlation", meansave)

# median
mediantitle = "Median absolute correlation of the factors through different number of factors\n(Adjusted by Fisher transformation)"
mediansave = "medianfishercorrperfactors.png"
linescatterplot("median", #2a9d8f, mediantitle, "Median of the correlation", mediansave)

# mode
modetitle = "Mode absolute correlation of the factors through different number of factors\n(Adjusted by Fisher transformation)"
modesave = "modefishercorrperfactors.png"
linescatterplot("mode", #e9c46a, modetitle , "Mode of the correlation", modesave )

