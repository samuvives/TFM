# the problem to solve is:
# the parameters we should manipulate easily are:
    # the number of factors
    # the views
    # the guiding variables

import os
import sys
import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import mofaflex as mfl
from mudata import MuData
from sklearn.preprocessing import StandardScaler
import anndata
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# =======================================================
# PARAMETERS
# =======================================================
print("Loading parameters")
# parameters es el modulo recibido como parametro que contiene los parametros del modelo
parameters = importlib.import_module(sys.argv[1])

INPUT_DIR = parameters.INPUT_DIR
WORKING_DIR = parameters.WORKING_DIR
N_FACTORS_LIST = parameters.N_FACTORS_LIST

VIEWS_CONFIG = parameters.VIEWS_CONFIG

GV_FILE = f"{INPUT_DIR}/matrizGV4_mapeado.csv"
BASE_OUTPUT = os.path.join(WORKING_DIR, "MOFAFLEX_FINAL_ANALYSIS")
os.makedirs(BASE_OUTPUT, exist_ok=True)

## funciones

# =======================================================
# 2. CARGA DE DATOS
# =======================================================
print("Loading views")
adata_dict = {}
for view, cfg in VIEWS_CONFIG.items():
    df = pd.read_csv(cfg["path"], sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if cfg["scale"]:
        df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)
    adata = anndata.AnnData(X=df.values.astype(np.float32), obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    adata_dict[view] = adata

# remove the last df used to liberate memory
del df 

mdata = MuData(adata_dict)
gv_df = pd.read_csv(GV_FILE, index_col=0)
gv_df.index = gv_df.index.astype(str)
gv_names = gv_df.columns.tolist()

for view in mdata.mod.keys():
    mdata.mod[view].obs = mdata.mod[view].obs.join(gv_df)


archivo_config = sys.argv[1]

# =======================================================
# 3 LOOP PRINCIPAL
# =======================================================
for k in N_FACTORS_LIST: 
    print(f"\n === K={k} ===")
    outdir = os.path.join(BASE_OUTPUT, f"K{k}")
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    model = mfl.MOFAFLEX(
        mdata,
        mfl.ModelOptions(
            likelihoods={v: VIEWS_CONFIG[v]["likelihood"] for v in VIEWS_CONFIG},
            guiding_vars_likelihoods={gv: "Bernoulli" for gv in gv_names},
            n_factors=k),
        mfl.TrainingOptions(seed=42),
        mfl.DataOptions(guiding_vars_obs_keys=gv_names))

    # ===================================================
    # ELBO
    # ===================================================
    elbofile = mfl.pl.training_curve(model)
    elbofile.save(os.path.join(outdir, f"training_curve_elbo_K{k}.png"), dpi=300)
    plt.close()

    # ===================================================
    # VARIANCE EXPLAINED PER VIEW
    # ===================================================
    # if total True
    # groups are columns and views are rows, and I only have one group
    varianceperview = model.get_r2(total=True)
    varianceperview.to_csv(os.path.join(outdir, f"varianceperview_K{k}.csv"))

    # ===================================================
    # VARIANCE EXPLAINED PER FACTOR
    # ===================================================
    # variance per factor method
    # total = False returns a dict of dataframe
    # each key is a group
    r2_dict = model.get_r2(total=False)
    nombre_grupo = list(r2_dict.keys())[0]
    df_r2 = r2_dict[nombre_grupo]
    # factors as rows, views as columns
    df_r2.to_csv(os.path.join(outdir, f"varianceperfactor_K{k}.csv"))

    # barplot stacked
    df_r2['Total_Var'] = df_r2.sum(axis=1)

    # 2. Ordenar de mayor a menor y eliminar la columna temporal para el gráfico
    df_r2 = df_r2.sort_values(by='Total_Var', ascending=False).drop(columns=['Total_Var'])
    df_r2_nonguided = df_r2[df_r2.index.str.startswith("Factor")]

    max_height = df_r2_nonguided.sum(axis=1).max()
    threshold = max_height * 0.01

    threshold2 = max_height * 0.02

    ax = df_r2.plot(kind='bar',
                  stacked=True,
                  figsize=(12, 7),
                  colormap="tab20",
                  edgecolor='white',
                  linewidth=0.5)

    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1, label=f'Threshold (1%)')
    ax.axhline(y=threshold2, color='green', linestyle='--', linewidth=1, label=f'Threshold (2%)')
    plt.title(f"Factor composition per view. K: {k}", fontsize=16, pad=20)
    plt.ylabel("Activity", fontsize=12)
    plt.xlabel("Factors", fontsize=12)
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")

    plt.legend(title="Omic views",
               bbox_to_anchor=(1.05, 1),
               loc='upper left',
               fontsize=10,
               frameon=False)

    plt.tight_layout()
    save_path = os.path.join(outdir, f"factor_composition_stacked_K{k}.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    # variance per factor heatmap
    varexpfile = mfl.pl.variance_explained(model)
    varexpfile.save(os.path.join(outdir, f"variance_explained_K{k}.png"), dpi=300)
    plt.close()

    # ===================================================
    # FULL MATRIX Z EXPORT
    # ===================================================
    # Matrix of samples x factors
    Z_dict = model.get_factors(ordered=True)
    group_name = list(Z_dict.keys())[0]
    df_Z_final = Z_dict[group_name]

    df_Z_final.to_csv(os.path.join(outdir, f"complete_factors_Z_K{k}.csv"))
    print("Full matrix Z exported")

    # ===================================================
    # FULL MATRICES W EXPORT PER VIEW
    # ===================================================
    # Matrices of features x factors
    weights = model.get_weights(ordered=True)

    weights_outdir = os.path.join(outdir, "complete_weights")
    os.makedirs(weights_outdir, exist_ok=True)

    for view, df_W in weights.items():
        # Transponemos (.T) para que: Filas = Genes/Variables, Columnas = Factores
        df_W_final = df_W.T

        file_name = f"complete_weights_{view}_K{k}.csv"
        df_W_final.to_csv(os.path.join(weights_outdir, file_name))

    print("Full matrices W exported")
    print(f"W.csv files saved in: {outdir}")

    print(f"Finalizado K={k}\n")

print("Program finished")
