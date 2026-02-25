# r2 is a measure of explained variance
# It allows to distinguish between general and specific factors

# MOFAFLEX.get_r2 total=False returns a dict of dataframes 
# variance explained for each view (columns) and factors (rows)

# MOFAFLEX.get_r2 total=True returns 
# each dataframe a group
# groups are columns and views are rows, and I only have one group
# where we have variance explained for each view (columns) and factors (rows)

# we need to do a heatmap and show the total variance too
# the threshold should be from the total
# it is not so simple as to obtain the total
# 2% of the variance in at least one of the views, makes the factor relevant
# you apply the threshold to all the views in a particular factor
# not to the sum of the variance of the views of that particular factor

#-----------------------------------------
# IMPORTS
#-----------------------------------------
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
pondsumlist_g = []
pondsumlist_ng = []

OUTPUTDIR = os.path.join(PATH, "postanalysis/r2analysis")
os.makedirs(OUTPUTDIR, exist_ok=True)

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

#-----------------------------------------
# FUNCTIONS
#-----------------------------------------
def stackedfactorsactivity(r2perfactor, KDIR, OUTPUTPATH, OUTPUTPATHWITHTHR):
    r2perfactor = r2perfactor.copy()

    # barplot stacked
    r2perfactor['Total_Var'] = r2perfactor.sum(axis=1)

    # 2. Ordenar de mayor a menor y eliminar la columna temporal para el gráfico
    r2perfactor = r2perfactor.sort_values(by='Total_Var', ascending=False).drop(columns=['Total_Var'])
    r2perfactor_nonguided = r2perfactor[r2perfactor.index.str.startswith("Factor")]

    max_height = r2perfactor_nonguided.sum(axis=1).max()
    threshold = max_height * 0.01
    threshold2 = max_height * 0.02

    ax = r2perfactor.plot(kind='bar',
                  stacked=True,
                  figsize=(12, 7),
                  colormap="tab20",
                  edgecolor='white',
                  linewidth=0.5)

    plt.title(f"Factor composition per view. {KDIR}", fontsize=16, pad=20)
    plt.ylabel("Activity", fontsize=12)
    plt.xlabel("Factors", fontsize=12)
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")
    plt.legend(title="Omic views",
               bbox_to_anchor=(1.05, 1),
               loc='upper left',
               fontsize=10,
               frameon=False)
    plt.tight_layout()
    plt.savefig(OUTPUTPATH, dpi=300, bbox_inches='tight')

    # graph with threshold
    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1, label=f'Threshold (1%)')
    ax.axhline(y=threshold2, color='green', linestyle='--', linewidth=1, label=f'Threshold (2%)')
    plt.savefig(OUTPUTPATHWITHTHR, dpi=300, bbox_inches='tight')

    plt.close()


# main graph
def maingraph(r2perfactor):
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
        r2perfactorpondsum[view] = r2perfactorpondsum[view] * totalfeaturesperview[view]

    # sum data and divide by total features
    r2perfactorpondsum["ponderedsumperfactor"] = r2perfactorpondsum.sum(axis=1) / totalfeatures
    r2perfactorpondsum = r2perfactorpondsum[["ponderedsumperfactor"]]
    return r2perfactorpondsum


def heatmap(viewsorder_numerico, viewsorderdf, mycolormap, r2perfactor, r2perfactorpondsum, r2perview, pathheatmap, KDIR):
    fig, axes = plt.subplots(2, 3,
        figsize=(20, 15),
        gridspec_kw = {
            "width_ratios": [8, 13, 1],
            "height_ratios": [r2perfactor.shape[0], 1]
            }
    )

    # upper left figure
    # columns = views ordered by relevance
    # rows = factors
    sns.heatmap(viewsorder_numerico,
                cmap=mycolormap,
                annot=viewsorderdf.values,
                fmt="",
                linewidths=0.5,
                linecolor='white',
                cbar=False,
                annot_kws = {"size": 8},
                ax=axes[0, 0])

    axes[0, 0].xaxis.set_ticks_position('top')
    axes[0, 0].xaxis.set_label_position('top')
    axes[0, 0].set_title("Views order by relevance in the factor")

    # upper mid figure
    # columns = views
    # rows = factors
    sns.heatmap(r2perfactor,
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

    # upper right figure
    # rows = relevance of each factor
    sns.heatmap(r2perfactorpondsum,
                cmap="Greens",
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

    axes[0, 2].set_title("Pondered sum")

    # bottom figure
    # columns = views
    sns.heatmap(r2perview,
                cmap="Blues",
                annot=True,
                linewidths=0.5, linecolor="white",
                ax=axes[1, 1])

    axes[1, 1].set_ylabel("")
    axes[1, 1].set_yticks([])
    axes[1, 1].set_yticklabels([])

    axes[1, 1].set_xlabel("")
    axes[1, 1].set_xticks([])
    axes[1, 1].set_xticklabels([])

    axes[1, 1].set_title("Views relevance in the model", y=-0.5)

    # eliminate lateral bottom figures
    fig.delaxes(axes[1, 0])
    fig.delaxes(axes[1, 2])

    fig.suptitle(f"Analysis of $R^2$: {KDIR}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(pathheatmap)
    plt.close()

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
# MAIN LOOP (GRAPHS PER RUN)
#-----------------------------------------
for KDIR in LISTDIRS:
    #-----------------------------------------
    # paths
    #-----------------------------------------
    path_KDIR = os.path.join(INPUTDIRPATH, KDIR)

    path_r2perfactor = os.path.join(path_KDIR, "varianza_explicada_por_factor.csv")
    r2perfactor = pd.read_csv(path_r2perfactor, index_col=0)

    path_r2perview = os.path.join(path_KDIR, "varianceperview.csv")

    #-----------------------------------------
    # stacked barplot
    #-----------------------------------------
    r2perfactor["Mean"] = r2perfactor.mean(axis=1)

    OUTPUTPATH = os.path.join(OUTPUTDIR, f"factor_composition_stacked_{KDIR}.png")
    OUTPUTPATHWITHTHR = os.path.join(OUTPUTDIR, f"factor_composition_stacked_thr_{KDIR}.png")
    stackedfactorsactivity(r2perfactor, KDIR, OUTPUTPATH, OUTPUTPATHWITHTHR)
    r2perfactor = r2perfactor.drop(columns=["Mean"])

    #-----------------------------------------
    # heatmap
    #-----------------------------------------
    # data for main heatmap
    # reduce the names
    r2perfactor = maingraph(r2perfactor)

    # data for bottom heatmap
    r2perview = getr2perview(path_r2perview, dictviewscorrected)

    # data for left heatmap
    viewsorder_numerico, viewsorderdf, mycolormap = orderviews(r2perfactor, colors_dict)

    # data for right heatmap
    r2perfactorpondsum = ponderedsum(r2perfactor, totalfeaturesperview, totalfeatures)

    # create the whole heatmap
    pathheatmap = os.path.join(OUTPUTDIR, f"r2_heatmap_{KDIR}.png")
    heatmap(viewsorder_numerico, viewsorderdf, mycolormap, r2perfactor, r2perfactorpondsum, r2perview, pathheatmap, KDIR)

    #-----------------------------------------
    # metrics to select the number of factors
    #-----------------------------------------
    # factors must explain at least 2% of the variance in at least one of the views.
    # is the r_2 a percentage?

    # create a column of the validity of the factors
    r2perfactor["FactorValidity"] = (r2perfactor > 0.02).any(axis=1)
    r2validityperfactor = r2perfactor[["FactorValidity"]]

    # r2perfactorpondsum has as index, the factors
    # "ponderedsumperfactor" as the single column
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

#-----------------------------------------
# GRAPHS COMPARING RUNS
#-----------------------------------------
# add non guided dataframes on top of each other
allpondsumdf_ng = pd.concat(pondsumlist_ng, axis=0)

# add guided dataframes on top of each other
allpondsumdf_g = pd.concat(pondsumlist_g, axis=0)

def r2_scatter_acrossruns(allpondsumdf_ng, allpondsumdf_g):
    """
    allpondsumdf_ng has 5 columns:
    Run, RunNumber, FactorNumber,
    ponderedsumperfactor, FactorValidity

    allpondsumdf_g has 5 columns:
    Run, RunNumber, FactorName,
    ponderedsumperfactor, FactorValidity
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
    pointcolors_ng = [mapcolors[c] for c in allpondsumdf_ng["FactorValidity"]]
    existing_runs = sorted(allpondsumdf_ng["RunNumber"].unique())

    axes[0].scatter(allpondsumdf_ng["RunNumber"], allpondsumdf_ng["FactorNumber"],
        c=pointcolors_ng,
        s=np.array(allpondsumdf_ng["ponderedsumperfactor"])*1000)

    axes[0].set_xticks(existing_runs)
    axes[0].set_xlabel("Run number")
    axes[0].set_ylabel("Factor number")
    axes[0].set_title("Variance explained in non-guided factors across the runs")

    # scatter for guided factors
    gf_forscatter = allpondsumdf_g["FactorName"].unique()
    map_gf = {cat: i for i, cat in enumerate(gf_forscatter)}
    allpondsumdf_g['y_num'] = allpondsumdf_g["FactorName"].map(map_gf)

    pointcolors_g = [mapcolors[c] for c in allpondsumdf_g["FactorValidity"]]

    axes[1].scatter(allpondsumdf_g["RunNumber"], allpondsumdf_g["y_num"],
        c=pointcolors_g,
        s=np.array(allpondsumdf_g["ponderedsumperfactor"])*1000)

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
    Line2D([0], [0], marker='o', color='w', label='Size of points = $r^2$',
           markerfacecolor='gray', markersize=10)
    ]
    fig.legend(handles=legend_elements,
           loc='lower left',
           bbox_to_anchor=(0.05, 0.02), # (X, Y) desde 0 a 1
           ncol=1,                       # ncol=1 para que los elementos vayan uno sobre otro
           frameon=True)

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig(os.path.join(OUTPUTDIR, "r2_scatter_acrossruns.png"))

r2_scatter_acrossruns(allpondsumdf_ng, allpondsumdf_g)
# revisar que esté todo cuadrado, que no se hayan dado la vuelta los datos
