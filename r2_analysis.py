# r2 is a measure of explained variance
# It allows to distinguish between general and specific factors

# MOFAFLEX.get_r2 total=False returns a dict of dataframes 
# variance explained for each view (columns) and factors (rows)

# MOFAFLEX.get_r2 total=True returns 
# each dataframe a group
# groups are columns and views are rows, and I only have one group
# where we have variance explained for each view (columns) and factors (rows)

# the threshold should be from the total
# it is not so simple as to obtain the total
# 2% of the variance in at least one of the views, makes the factor relevant
# you apply the threshold to all the views in a particular factor
# not to the sum of the variance of the views of that particular factor
# I need to divide the main heatmap in two parts separating the guided and nonguided factors

#-----------------------------------------
# IMPORTS
#-----------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D

import seaborn as sns

#-----------------------------------------
# PARAMETERS
#-----------------------------------------
approach = "OBTAININGELBO"

PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/{approach}"

INPUTDIRPATH = os.path.join(PATH, "MOFAFLEX_FINAL_ANALYSIS")
LISTDIRS = [f for f in os.listdir(INPUTDIRPATH) if os.path.isdir(os.path.join(INPUTDIRPATH, f))]
pathtotalfeaturesperview = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/featuresandsparsity/tableuniquevaluessparsity.csv"
INTLISTDIRS = [int(f.replace("K", "")) for f in LISTDIRS]
maxrunnumber = max(INTLISTDIRS)
pondsumlist_ng = []
pondsumlist_g = []

OUTPUTDIR = os.path.join(PATH, "postanalysis/r2analysis")
os.makedirs(OUTPUTDIR, exist_ok=True)

# dict for reducing the names of the views
dictviewscorrected = {
    "Metabolomics": "MET",
    "Lipidomics": "LIP",
    "Microbiota": "MIC",
    "SV_DEL": "SV_DEL",
    "SV_DUP": "SV_DUP",
    "SV_INS": "SV_INS",
    "SV_INV": "SV_INV",
    "SV_TRA": "SV_TRA",
    "VC_11": "VC_11",
    "VC_12": "VC_12",
    "EXPRESSION": "EXP"
}

# dicts for the colors in the left heatmaps
colors_dict = {
    "Metabolomics": "#99582a",
    "Lipidomics": "#ffd23f",
    "Microbiota": "grey",
    "SV_DEL": "#3bceac",
    "SV_DUP": "#0496ff",
    "SV_INS": "#007ea7",
    "SV_INV": "#3da5d9",
    "SV_TRA": "#758bfd",
    "VC_11": "#c77dff",
    "VC_12": "#5a189a",
    "EXPRESSION": "#ee4266"
}

dicttotalr2 = {}

#-----------------------------------------
# FUNCTIONS
#-----------------------------------------

# main graph
def datamaingraph(r2perfactor):
    r2perfactor = r2perfactor.rename(columns=dictviewscorrected)
    return r2perfactor

# bottom graph: variance per view
def getr2perview(path_r2perview, dictviewscorrected):
    r2perview = pd.read_csv(path_r2perview, index_col=0)
    r2perview.index = list(dictviewscorrected.values())
    r2perview = r2perview.T
    return r2perview

# left graph: relevance of omics in each factor

def orderviews(r2perfactor, colors_dict):
    viewsorderdict = {}
    for factor in r2perfactor.index.tolist():
        viewsorder = r2perfactor.loc[factor].sort_values(ascending=False).index.tolist()
        viewsorderdict[factor] = viewsorder
    viewsorderdf = pd.DataFrame(viewsorderdict)
    # each column a factor, each row a position
    viewsorderdf = viewsorderdf.T
    viewsorderdf.columns = [i for i in range(1, viewsorderdf.shape[1] + 1)]
    # each column a position, each row a factor

    # reduce the names of colors_dict
    colors_dict = {dictviewscorrected[k]: v for k, v in colors_dict.items()}

    # stablish num matrix for seaborn to understand
    vistas_unicas = list(colors_dict.keys())
    vista_a_numero = {vista: i for i, vista in enumerate(vistas_unicas)}

    # Crear matriz numérica para los COLORES (no para mostrar)
    viewsorder_numerico = viewsorderdf.map(lambda x: vista_a_numero[x])

    # Crear lista de colores en el mismo orden
    colores_lista = [colors_dict[vista] for vista in vistas_unicas]
    mycolormap = ListedColormap(colores_lista)
    return viewsorder_numerico, viewsorderdf, mycolormap

# right graph: pondered sum
# inside the loop
def ponderedsum(r2perfactor, totalfeaturesperview, totalfeatures) -> pd.DataFrame:
    # make a copy of the data to create a table of pondered sum
    r2perfactorpondsum = r2perfactor.copy()

    # multiply data by total features per view
    for view in r2perfactorpondsum.columns:
        r2perfactorpondsum[view] = r2perfactorpondsum[view]

    # sum data and divide by total features
    r2perfactorpondsum["Sum"] = r2perfactorpondsum.sum(axis=1)
    r2perfactorpondsum = r2perfactorpondsum[["Sum"]]
    return r2perfactorpondsum


def heatmap(
    viewsorder_numerico_ng, viewsorder_numerico_g,
    viewsorderdf_ng, viewsorderdf_g,
    mycolormap_ng, mycolormap_g, 
    r2perfactor_ng, r2perfactor_g,
    r2perfactorpondsum_ng, r2perfactorpondsum_g,
    r2perview, pathheatmap, KDIR):

    fig, axes = plt.subplots(3, 3,
        figsize=(20, r2perfactor_ng.shape[0]),
        gridspec_kw = {
            "width_ratios": [2, 13, 1],
            "height_ratios": [r2perfactor_ng.shape[0], r2perfactor_g.shape[0], 1]
            }
    )

    # upper left figures
    # columns = views ordered by relevance
    # rows = factors
    sns.heatmap(viewsorder_numerico_ng,
                cmap=mycolormap_ng,
                linewidths=0.5,
                linecolor='white',
                cbar=False,
                annot_kws = {"size": 8},
                ax=axes[0, 0])

    axes[0, 0].xaxis.set_ticks_position('top')
    axes[0, 0].xaxis.set_label_position('top')
    axes[0, 0].set_title("Views order by relevance in the factor")
    # create legend of the colors, and do not put the annotations in the boxes
    patches = [mpatches.Patch(color=color, label=nombre) for nombre, color in colors_dict.items()]
    axes[0, 0].legend(
        handles=patches, 
        bbox_to_anchor=(-0.3, 0.5),
        loc="upper right",
        borderaxespad=0.,
        frameon=False)

    sns.heatmap(viewsorder_numerico_g,
                cmap=mycolormap_g,
                linewidths=0.5,
                linecolor='white',
                cbar=False,
                annot_kws = {"size": 8},
                ax=axes[1, 0])

    axes[1, 0].set_xlabel("")
    axes[1, 0].set_xticks([])
    axes[1, 0].set_xticklabels([])

    # upper mid figure
    # columns = views
    # rows = factors
    sns.heatmap(r2perfactor_ng,
                cmap="Reds",
                annot=True,
                linewidths=0.5,
                linecolor="white",
                cbar_kws={"shrink": 0.6},
                ax=axes[0, 1])

    axes[0, 1].xaxis.set_ticks_position("top")
    axes[0, 1].xaxis.set_label_position("top")

    axes[0, 1].set_ylabel("")
    axes[0, 1].set_yticks([])
    axes[0, 1].set_yticklabels([])

    axes[0, 1].set_title("$R^2$")

    sns.heatmap(r2perfactor_g,
                cmap="Greens",
                annot=True,
                linewidths=0.5,
                linecolor="white",
                cbar_kws={"shrink": 0.6},
                ax=axes[1, 1])

    axes[1, 1].set_xlabel("")
    axes[1, 1].set_xticks([])
    axes[1, 1].set_xticklabels([])

    axes[1, 1].set_ylabel("")
    axes[1, 1].set_yticks([])
    axes[1, 1].set_yticklabels([])

    # upper right figure
    # rows = relevance of each factor
    sns.heatmap(r2perfactorpondsum_ng,
                cmap="Reds",
                annot=True,
                linewidths=0.5, 
                linecolor="white",
                ax=axes[0, 2])

    axes[0, 2].set_xlabel("")
    axes[0, 2].set_xticks([])
    axes[0, 2].set_xticklabels([])

    axes[0, 2].set_ylabel("")
    axes[0, 2].set_yticks([])
    axes[0, 2].set_yticklabels([])

    axes[0, 2].set_title("Sum")

    sns.heatmap(r2perfactorpondsum_g,
                cmap="Greens",
                annot=True,
                linewidths=0.5, 
                linecolor="white",
                ax=axes[1, 2])

    axes[1, 2].set_xlabel("")
    axes[1, 2].set_xticks([])
    axes[1, 2].set_xticklabels([])

    axes[1, 2].set_ylabel("")
    axes[1, 2].set_yticks([])
    axes[1, 2].set_yticklabels([])

    # bottom figure
    # columns = views
    sns.heatmap(r2perview,
                cmap="Blues",
                annot=True,
                linewidths=0.5,
                linecolor="white",
                ax=axes[2, 1])

    axes[2, 1].set_ylabel("")
    axes[2, 1].set_yticks([])
    axes[2, 1].set_yticklabels([])

    axes[2, 1].set_xlabel("")
    axes[2, 1].set_xticks([])
    axes[2, 1].set_xticklabels([])

    axes[2, 1].set_title("Views relevance in the model", y=-0.5)

    # eliminate lateral bottom figures
    fig.delaxes(axes[2, 0])
    fig.delaxes(axes[2, 2])

    fig.suptitle(f"Analysis of $R^2$: {KDIR}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(pathheatmap)
    plt.close()

def r2_scatter_acrossruns(allpondsumdf_ng, allpondsumdf_g):
    """
    Receive dfs with all the factors from all the runs
    allpondsumdf_ng has 5 columns:
    Run, RunNumber, FactorNumber,
    Sum, FactorValidity

    allpondsumdf_g has 5 columns:
    Run, RunNumber, FactorName,
    Sum, FactorValidity
    """
    # color map for validity of factors
    mapcolors = {False: "red", True: "green"}

    # figure with 2 rows, 1 column
    fig, axes = plt.subplots(2, 1, figsize=(15, 15),
        gridspec_kw = {
            "height_ratios": [10, 4]
            }
    )

    # scatter for non guided factors
    # make the cmap
    pointcolors_ng = [mapcolors[c] for c in allpondsumdf_ng["FactorValidity"]]
    existing_runs = sorted(allpondsumdf_ng["RunNumber"].unique())

    axes[0].scatter(allpondsumdf_ng["RunNumber"], allpondsumdf_ng["FactorNumber"],
        c=pointcolors_ng,
        s=np.array(allpondsumdf_ng["Sum"])*1000)

    axes[0].set_xticks(existing_runs)
    axes[0].set_xlabel("Run number")
    axes[0].set_ylabel("Factor number")
    axes[0].set_title("Variance explained in non-guided factors across the runs")

    # scatter for guided factors
    # map to points to visualize
    gf_forscatter = allpondsumdf_g["FactorName"].unique()
    map_gf = {cat: i for i, cat in enumerate(gf_forscatter)}
    allpondsumdf_g['y_num'] = allpondsumdf_g["FactorName"].map(map_gf)

    # make the cmap
    pointcolors_g = [mapcolors[c] for c in allpondsumdf_g["FactorValidity"]]

    axes[1].scatter(allpondsumdf_g["RunNumber"], allpondsumdf_g["y_num"],
        c=pointcolors_g,
        s=np.array(allpondsumdf_g["Sum"])*1000)

    axes[1].set_yticks(list(map_gf.values()))
    axes[1].set_yticklabels(list(map_gf.keys()))

    axes[1].set_xticks(existing_runs)
    axes[1].set_xlabel("Run number")
    axes[1].set_ylim(-1, 2)
    axes[1].set_ylabel("Factor number")
    axes[1].set_title("Variance explained in guided factors across the runs")

    fig.suptitle("Variance explained (pondered sum of the views $r^2$) \n of the factors across the runs", fontsize=13, fontweight="bold")

    legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Valid (>2% in at least 1 view)',
           markerfacecolor='green', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Not Valid',
           markerfacecolor='red', markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Size of points = $Sum r^2$',
           markerfacecolor='gray', markersize=10)
    ]
    fig.legend(handles=legend_elements,
           loc='lower left',
           bbox_to_anchor=(0.05, 0.02), # (X, Y) desde 0 a 1
           ncol=1,                       # ncol=1 para que los elementos vayan uno sobre otro
           frameon=True)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig(os.path.join(OUTPUTDIR, "r2_scatter_acrossruns.png"))

# def r2_heatmap_acrossruns(allpondsumdf_ng, allpondsumdf_g):


#-----------------------------------------
# right graph: total input features
#-----------------------------------------
# outside the loop
totalfeaturesperview = pd.read_csv(pathtotalfeaturesperview)
totalfeatures = totalfeaturesperview["Total number of features"].sum()
# select the two columns, convert to a dict
totalfeaturesperview = dict(zip(totalfeaturesperview["View"], totalfeaturesperview["Total number of features"]))
totalfeaturesperview = {dictviewscorrected[k]: v for k, v in totalfeaturesperview.items()}

#-----------------------------------------
# lists of series for validity heatmap
#-----------------------------------------
list_factorseries_ng_sum = []
list_factorseries_ng_valid = []
list_factorseries_g_sum = []
list_factorseries_g_valid = []
#-----------------------------------------
# MAIN LOOP (GRAPHS PER RUN)
#-----------------------------------------
for KDIR in LISTDIRS:
    #-----------------------------------------
    # paths
    #-----------------------------------------
    path_KDIR = os.path.join(INPUTDIRPATH, KDIR)

    path_r2perfactor = os.path.join(path_KDIR, "varianza_explicada_por_factor.csv")
    r2perfactor = pd.read_csv(path_r2perfactor, index_col=0)
    # rows factors, views columns

    path_r2perview = os.path.join(path_KDIR, "varianceperview.csv")
    # rows views, groups columns
    # only column: group_1

    #-----------------------------------------
    # heatmap
    #-----------------------------------------
    # data for 2 center heatmap
    # reduce the names
    r2perfactor = datamaingraph(r2perfactor)
    print("-" * 20)
    print("r2perfactor: ")
    print(r2perfactor)
    print("-" * 20)
    # divide between guided and non-guided factors
    r2perfactor_ng = r2perfactor.loc[[row for row in r2perfactor.index if row.startswith("Factor")]].copy()
    print("r2perfactor_ng:")
    print(r2perfactor_ng)
    print("-" * 20)
    r2perfactor_g = r2perfactor.loc[[row for row in r2perfactor.index if not row.startswith("Factor")]].copy()
    print("r2perfactor_g:")
    print(r2perfactor_g)

    # data for 2 left heatmap
    viewsorder_numerico_ng, viewsorderdf_ng, mycolormap_ng = orderviews(r2perfactor_ng, colors_dict)
    viewsorder_numerico_g, viewsorderdf_g, mycolormap_g = orderviews(r2perfactor_g, colors_dict)

    # data for 2 right heatmap
    r2perfactorpondsum = ponderedsum(r2perfactor, totalfeaturesperview, totalfeatures)
    # only one column "Sum"

    r2perfactorpondsum_ng = ponderedsum(r2perfactor_ng, totalfeaturesperview, totalfeatures)

    print("r2perfactorpondsum_ng:")
    print(r2perfactorpondsum_ng)
    r2perfactorpondsum_g = ponderedsum(r2perfactor_g, totalfeaturesperview, totalfeatures)
    print("r2perfactorpondsum_g:")
    print(r2perfactorpondsum_g)

    # data for bottom heatmap
    r2perview = getr2perview(path_r2perview, dictviewscorrected)
    # one row, each column a view

    # total r2
    total_r2 = r2perview.sum(axis=1)
    print("-" * 20)
    print("Total_r2: ")
    print(total_r2)
    dicttotalr2[KDIR] = total_r2

    # create the whole heatmap
    pathheatmap = os.path.join(OUTPUTDIR, f"r2_heatmap_{KDIR}.png")
    heatmap(viewsorder_numerico_ng, viewsorder_numerico_g,
        viewsorderdf_ng, viewsorderdf_g,
        mycolormap_ng, mycolormap_g, 
        r2perfactor_ng, r2perfactor_g,
        r2perfactorpondsum_ng, r2perfactorpondsum_g,
        r2perview, pathheatmap, KDIR)

    #-----------------------------------------
    # metrics to select the number of factors
    #-----------------------------------------
    # factors must explain at least 2% of the variance in at least one of the views.
    # is the r_2 a percentage?

    # create a column of the validity of the factors
    r2perfactor["FactorValidity"] = (r2perfactor > 0.02).any(axis=1)
    r2validityperfactor = r2perfactor[["FactorValidity"]]

    # r2perfactorpondsum has as index, the factors
    # "Sum" as the single column
    # add the validity column
    r2perfactorpondsum = pd.concat([r2perfactorpondsum, r2validityperfactor], axis=1)

    # create a column of the name of the run
    r2perfactorpondsum["Run"] = [KDIR] * r2perfactorpondsum.shape[0]
    r2perfactorpondsum["RunNumber"] = r2perfactorpondsum["Run"].str.replace("K", "").astype(int)

    # divide the df between guided and nonguided (divide by index)
    r2perfactorpondsum_ng = r2perfactorpondsum[r2perfactorpondsum.index.str.startswith("Factor")].copy()
    r2perfactorpondsum_g = r2perfactorpondsum[~r2perfactorpondsum.index.str.startswith("Factor")].copy()

    # create FactorNumber column for non guided factors df
    r2perfactorpondsum_ng["FactorNumber"] = r2perfactorpondsum_ng.index
    r2perfactorpondsum_ng["FactorNumber"] = r2perfactorpondsum_ng["FactorNumber"].str.replace("Factor ", "").astype(int)

    # create FactorName column for guided factors df
    r2perfactorpondsum_g["FactorName"] = r2perfactorpondsum_g.index

    # add the df to a list of dfs
    pondsumlist_ng.append(r2perfactorpondsum_ng)
    pondsumlist_g.append(r2perfactorpondsum_g)
    
    #--------
    # add the df to a list of dfs for validity heatmaps
    # fill the df with na if needed
    # NON GUIDED FACTORS
    r2perfactorvalidity_ng = r2perfactorpondsum_ng.copy()
    r2perfactorvalidity_ng = r2perfactorvalidity_ng.set_index("FactorNumber")
    r2perfactorvalidity_ng = r2perfactorvalidity_ng.reindex(range(1, maxrunnumber + 1))

    # sum heatmap
    factorseries_ng_sum = r2perfactorvalidity_ng["Sum"]
    factorseries_ng_sum.name = KDIR
    list_factorseries_ng_sum.append(factorseries_ng_sum)

    # validity heatmap
    factorseries_ng_valid = r2perfactorvalidity_ng["FactorValidity"]
    factorseries_ng_valid.name = KDIR
    list_factorseries_ng_valid.append(factorseries_ng_valid)

    # GUIDED FACTORS
    r2perfactorvalidity_g = r2perfactorpondsum_g.copy()
    r2perfactorvalidity_g = r2perfactorvalidity_g.set_index("FactorName")

    # sum heatmap
    factorseries_g_sum = r2perfactorvalidity_g["Sum"]
    factorseries_g_sum.name = KDIR
    list_factorseries_g_sum.append(factorseries_g_sum)

    # validity heatmap
    factorseries_g_valid = r2perfactorvalidity_g["FactorValidity"]
    factorseries_g_valid.name = KDIR
    list_factorseries_g_valid.append(factorseries_g_valid)

#-----------------------------------------
# GRAPHS COMPARING RUNS
#-----------------------------------------
# scatter

# add non guided dataframes on top of each other
allpondsumdf_ng = pd.concat(pondsumlist_ng, axis=0)

# add guided dataframes on top of each other
allpondsumdf_g = pd.concat(pondsumlist_g, axis=0)

r2_scatter_acrossruns(allpondsumdf_ng, allpondsumdf_g)

# df for the heatmaps
# non guided
dfvalidityheatmap_ng_s = pd.concat(list_factorseries_ng_sum, axis=1)
dfvalidityheatmap_ng_s = dfvalidityheatmap_ng_s.sort_index(ascending=False)
dfvalidityheatmap_ng_s.columns = [int(f.replace("K", "")) for f in dfvalidityheatmap_ng_s.columns]
dfvalidityheatmap_ng_s = dfvalidityheatmap_ng_s.sort_index(axis=1)

dfvalidityheatmap_ng_val = pd.concat(list_factorseries_ng_valid, axis=1)
dfvalidityheatmap_ng_val = dfvalidityheatmap_ng_val.fillna(0).astype(int)
dfvalidityheatmap_ng_val = dfvalidityheatmap_ng_val.sort_index(ascending=False)
dfvalidityheatmap_ng_val.columns = [int(f.replace("K", "")) for f in dfvalidityheatmap_ng_val.columns]
dfvalidityheatmap_ng_val = dfvalidityheatmap_ng_val.sort_index(axis=1)

# guided
dfvalidityheatmap_g_s = pd.concat(list_factorseries_g_sum, axis=1)
dfvalidityheatmap_g_s.columns = [int(f.replace("K", "")) for f in dfvalidityheatmap_g_s.columns]
dfvalidityheatmap_g_s = dfvalidityheatmap_g_s.sort_index(axis=1)

dfvalidityheatmap_g_val = pd.concat(list_factorseries_g_valid, axis=1)
dfvalidityheatmap_g_val = dfvalidityheatmap_g_val.fillna(0).astype(int)
dfvalidityheatmap_g_val.columns = [int(f.replace("K", "")) for f in dfvalidityheatmap_g_val.columns]
dfvalidityheatmap_g_val = dfvalidityheatmap_g_val.sort_index(axis=1)

# 2 heatmaps overlapping of ng factors
#"Run", "FactorValidity", "Sum"
binarycmap = ListedColormap(["#FFFFFF00", "#00FF00"])
fig, axes = plt.subplots(2, 1, figsize = (16, 16),
            gridspec_kw = {
                "height_ratios": [10, 4]
            })
sns.heatmap(data = dfvalidityheatmap_ng_s,
            cmap="Reds",
            annot=True,
            fmt=".3f",
            linewidths=0.5, 
            linecolor="white",
            ax=axes[0])

sns.heatmap(data = dfvalidityheatmap_ng_val,
            mask=(dfvalidityheatmap_ng_val == 0),
            cmap=binarycmap,
            alpha=0.8,
            cbar=False,
            annot=False,
            linewidths=0.5, 
            linecolor="white",
            ax=axes[0],
            zorder=10)

axes[0].set_xlabel("")
axes[0].set_xticks([])
axes[0].set_xticklabels([])
axes[0].set_title("Non Guided Factors")

# 2 heatmaps overlapping of g factors
sns.heatmap(data = dfvalidityheatmap_g_s,
            cmap="Reds",
            annot=True,
            fmt=".3f",
            linewidths=0.5, 
            linecolor="white",
            ax=axes[1])

sns.heatmap(data = dfvalidityheatmap_g_val,
            mask=(dfvalidityheatmap_g_val == 0),
            cmap=binarycmap,
            alpha=0.8,
            cbar=False,
            annot=False,
            linewidths=0.5, 
            linecolor="white",
            ax=axes[1],
            zorder=10)

axes[1].set_title("Guided Factors")
fig.suptitle(f"Analysis of validity of factors \n ($R^2$) through runs", fontsize=13, fontweight="bold")
plt.savefig(os.path.join(OUTPUTDIR, "validityheatmap.png"))
plt.close()

# total r2
dftotalr2 = pd.DataFrame.from_dict(dicttotalr2, orient="index").reset_index()
print("df total r2: ")
print(dftotalr2)
dftotalr2.columns = ["Run", "Total_r2"]
dftotalr2["RunNumber"] = dftotalr2["Run"].str.replace("K", "").astype(int)
dftotalr2 = dftotalr2.sort_values(by="RunNumber")
fig, ax = plt.subplots()
ax.scatter(dftotalr2["RunNumber"], dftotalr2["Total_r2"], color= "#2a9d8f")
ax.plot(dftotalr2["RunNumber"], dftotalr2["Total_r2"], color= "#2a9d8f")
ax.set_xlabel("Run with k factors")
ax.set_ylabel("Total r2")
ax.set_xticks(range(len(dftotalr2["RunNumber"]))
ax.set_xticklabels(dftotalr2["RunNumber"])
ax.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUTDIR, "totalr2lineplot.png"))
plt.close()
