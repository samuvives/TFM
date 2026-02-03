import os
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

## funciones
def zfactors(file, outdir, gv_names):
    """Funcion para saber los factores inutiles"""
    fdata = pd.read_csv(file, index_col=0)
    # Filtrar solo los que empiezan por "Factor" para evitar las GVs si estuvieran ahí
    fdata = fdata[fdata.index.str.contains("Factor")]
    
    plt.figure(figsize=(10, 5))
    sns.barplot(x=fdata.index, y=fdata.iloc[:, 0], palette="viridis")
    plt.ylabel("Activity")
    plt.xlabel("Factors")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "activityNGfactors.png"), dpi=300)

    max_valor = fdata.iloc[:, 0].max()
    thr = 0.01 * max_valor
    mask = fdata.iloc[:, 0] < thr
    inactive_list = fdata.index[mask].tolist()
    if len(inactive_list) == 0:
        print("Todos los factores superan el umbral de actividad")
    else:
        print(f"Los factores inactivos son: {', '.join(inactive_list)}")


def factoractivityperview(weights):
    # Calculamos la actividad de cada vista por factor (Suma de cuadrados de pesos)
    # Esto indica cuánto contribuye cada ómica a cada factor
    # elevas al cuadrado todas las celdas y sumas a traves de las columnas, quedandote solo con una columna
    # cada fila de esa columna es un factor
    # para cada view tienes un panda series con el valor de cada factor
    # eso lo pones en un diccionario, que conviertes en dataframe
    # cada columna es una view, cada fila un factor
    view_activity_dict = {view: (W**2).sum(axis=1) for view, W in weights.items()}
    df_activity_all = pd.DataFrame(view_activity_dict)
    
    # Ordenar los factores por actividad total (de mayor a menor)
    # sumas todas las columnas y te quedas con un valor global por factor
    total_activity = df_activity_all.sum(axis=1).sort_values(ascending=False)
    df_plot = df_activity_all.loc[total_activity.index]

    return df_plot


# =======================================================
# 1. PARÁMETROS
# =======================================================
print("[1] Inicializando parámetros")
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal"
N_FACTORS_LIST = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

VIEWS_CONFIG = {
    "Metabolomics": {"path": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Lipidomics": {"path": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Microbiota": {"path": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv", "likelihood": "Normal", "scale": True},
    "SV_DEL": {"path": f"{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_DUP": {"path": f"{INPUT_DIR}/SV_DUP_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_INS": {"path": f"{INPUT_DIR}/SV_INS_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_INV": {"path": f"{INPUT_DIR}/SV_INV_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "SV_TRA": {"path": f"{INPUT_DIR}/SV_TRA_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
    "VC_11": {"path": f"{INPUT_DIR}/MPA_GT_1_1.tsv", "likelihood": "Bernoulli", "scale": False},
    "VC_12": {"path": f"{INPUT_DIR}/MPA_GT_1_2.tsv", "likelihood": "Bernoulli", "scale": False},
    "EXPRESSION": {"path": f"{INPUT_DIR}/tpmexpression.tsv", "likelihood": "Normal", "scale": True}
}

GV_FILE = f"{INPUT_DIR}/matrizGV4_mapeado.csv"
BASE_OUTPUT = os.path.join(WORKING_DIR, "MOFAFLEX_FINAL_ANALYSIS")
os.makedirs(BASE_OUTPUT, exist_ok=True)

# =======================================================
# 2. CARGA DE DATOS
# =======================================================
print("[2] Cargando vistas")
adata_dict = {}
for view, cfg in VIEWS_CONFIG.items():
    df = pd.read_csv(cfg["path"], sep="\t", index_col=0)
    df.index = df.index.astype(str)
    if cfg["scale"]:
        df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)
    adata = anndata.AnnData(X=df.values.astype(np.float32), obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    adata_dict[view] = adata

mdata = MuData(adata_dict)
gv_df = pd.read_csv(GV_FILE, index_col=0)
gv_df.index = gv_df.index.astype(str)
gv_names = gv_df.columns.tolist()

for view in mdata.mod.keys():
    mdata.mod[view].obs = mdata.mod[view].obs.join(gv_df)

# =======================================================
# 3 LOOP PRINCIPAL
# =======================================================
for k in N_FACTORS_LIST: # meter aqui un generator?
    print(f"\n === K={k} ===")
    outdir = os.path.join(BASE_OUTPUT, f"K{k}")
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    auto_save_path = os.path.join(outdir, f"model_K{k}.pkl")

    model = mfl.MOFAFLEX(
        mdata,
        mfl.ModelOptions(
            likelihoods={v: VIEWS_CONFIG[v]["likelihood"] for v in VIEWS_CONFIG},
            guiding_vars_likelihoods={gv: "Bernoulli" for gv in gv_names},
            n_factors=k),
        mfl.TrainingOptions(seed=42, save_path=auto_save_path),
        mfl.DataOptions(guiding_vars_obs_keys=gv_names))

    # ===================================================
    # ELBO
    # ===================================================
    elbofile = mfl.pl.training_curve(model)
    elbofile.save(os.path.join(outdir, "training_curve_elbo.png"), dpi=300)

    # ===================================================
    # VARIANCE EXPLAINED
    # ===================================================
    varexpfile = mfl.pl.variance_explained(model)
    varexpfile.save(os.path.join(outdir, "variance_explained.png"), dpi=300)

    # ===================================================
    # FULL MATRIX Z EXPORT
    # ===================================================
    # Usamos el método nativo para recuperar la tabla de muestras x factores
    Z_dict = model.get_factors(ordered=True)
    group_name = list(Z_dict.keys())[0]
    df_Z_final = Z_dict[group_name]

    df_Z_final.to_csv(os.path.join(outdir, f"complete_factors_Z_K{k}.csv"))
    print("Full matrix Z exported")

    # ===================================================
    # FULL MATRICES W EXPORT PER VIEW
    # ===================================================
    # Usamos el método nativo para recuperar la tabla de variables x factores
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

    # ===================================================
    # 3.1 Z Activity
    # ===================================================
    Z = Z_dict[list(Z_dict.keys())[0]]
    factor_activity = (Z ** 2).mean(axis=0)

    zfactorsfile = os.path.join(outdir, "factor_activity.csv")
    factor_activity.to_csv(zfactorsfile)

    zfactors(zfactorsfile, outdir, gv_names)

    # ===================================================
    # 3.2 Factors activity using Weights & Stacked Bar
    # ===================================================
    from stackedmodified import stackedgraph
    
    # Calculamos la actividad de cada vista por factor (Suma de cuadrados de pesos)
    # Esto indica cuánto contribuye cada ómica a cada factor
    view_activity_dict = {view: (W**2).sum(axis=1) for view, W in weights.items()}
    df_activity_all = pd.DataFrame(view_activity_dict)
    
    # Ordenar los factores por actividad total (de mayor a menor)
    total_activity = df_activity_all.sum(axis=1).sort_values(ascending=False)
    df_plot = df_activity_all.loc[total_activity.index]

    stackedgraph(df_plot)

    # ===================================================
    # 3.14 CORRELACION CLINICA (FACTORES VS GV)
    # ===================================================
    print("[3.14] Calculando correlación clínica")

    from clinicalcorrelationmodified import clinicalcorrelationheatmap
    INPUT_BASE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/"
    INPUT_DIR = os.path.join(INPUT_BASE, k_number)
    GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"
    OUTPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis"

    # load factors
    factors_file = os.path.join(INPUT_DIR, f"complete_factors_Z_{k_number}.csv")
    Z = pd.read_csv(factors_file, index_col=0)

    print("Z head is:")
    print(Z.head(3))

    # load gv
    gv_df = pd.read_csv(GV_FILE, index_col=0)
    gv_names = gv_df.columns.tolist()

    print("gv file head is:")
    print(gv_df.head(3))
    clinicalcorrelationheatmap(Z, gv_names, k, OUTPUT_DIR)

    print(f"Finalizado K={k}\n")

print("Program finished")
