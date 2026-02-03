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


# =======================================================
# 1. PARÁMETROS
# =======================================================
print("[1] Inicializando parámetros")
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUTSUBSETS"
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/APPROACHSUBSETS"
N_FACTORS_LIST = [10, 12]

VIEWS_CONFIG = {
    "Metabolomics": {"path": f"{INPUT_DIR}/renamed_Metabolomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Lipidomics": {"path": f"{INPUT_DIR}/renamed_lipidomics_data_case.tsv", "likelihood": "Normal", "scale": True},
    "Microbiota": {"path": f"{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv", "likelihood": "Normal", "scale": True},
    "SV_DEL": {"path": f"{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv", "likelihood": "Bernoulli", "scale": False},
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
for k in N_FACTORS_LIST:
    print(f"\n === K={k} ===")
    outdir = os.path.join(BASE_OUTPUT, f"K{k}")
    figdir = os.path.join(outdir, "figures")
    os.makedirs(figdir, exist_ok=True)
    auto_save_path = os.path.join(outdir, f"model_K{k}.pkl")

    model = mfl.MOFAFLEX(mdata, mfl.ModelOptions(
            likelihoods={v: VIEWS_CONFIG[v]["likelihood"] for v in VIEWS_CONFIG},
            guiding_vars_likelihoods={gv: "Bernoulli" for gv in gv_names},
            n_factors=k),
        mfl.TrainingOptions(seed=42, save_path=auto_save_path),
        mfl.DataOptions(guiding_vars_obs_keys=gv_names))

    # ===================================================
    # ELBO
    # ===================================================
    # ELBO plot
    elbofile = mfl.pl.training_curve(model)
    elbofile.save(os.path.join(outdir, "training_curve_elbo.png"), dpi=300)

    # ELBO values (training epochs)
    elbo_values = model.training_loss
    print(elbo_values[:10])

    # ===================================================
    # VARIANCE EXPLAINED PER VIEW
    # ===================================================
    varianceperview = model.get_r2(total=True)
    varianceperview.to_csv(os.path.join(outdir,"varianceperview.csv"))
    print("head of the variance explained per view file:")
    print(head(varianceperview))

    # ===================================================
    # VARIANCE EXPLAINED PER FACTOR
    # ===================================================
    # variance per factor values
    r2_dict = model.get_r2(total=False)
    nombre_grupo = list(r2_dict.keys())[0]
    df_r2 = r2_dict[nombre_grupo]
    # factors as rows, views as columns
    print("head of the variance explained per factor file:")
    print(head(df_r2))
    df_r2.to_csv(os.path.join(outdir,"varianza_explicada_por_factor.csv"))

    # variance per factor plot
    varexpfile = mfl.pl.variance_explained(model)
    varexpfile.save(os.path.join(outdir, "variance_explained.png"), dpi=300)

    # ===================================================
    # TOP WEIGHTS
    # ===================================================
    topweightsfile = mfl.pl.top_weights(model)
    topweightsfile.save(os.path.join(outdir, "topweightsperfactor.png"), dpi=300)

    # ===================================================
    # 3.15 EXPORTACIÓN DE MATRICES COMPLETAS (Z y W)
    # ===================================================
    print(f"[3.15] Exportando matrices completas para análisis externo (K={k})")

    # 1. Obtener y guardar FACTORES (Z) completos
    # Usamos el método nativo para recuperar la tabla de muestras x factores
    Z_complete_dict = model.get_factors(ordered=True)
    group_name = list(Z_complete_dict.keys())[0]
    df_Z_final = Z_complete_dict[group_name]

    df_Z_final.to_csv(os.path.join(outdir, f"complete_factors_Z_K{k}.csv"))

    # 2. Obtener y guardar PESOS (W) completos por cada vista
    # Usamos el método nativo para recuperar la tabla de variables x factores
    W_complete_dict = model.get_weights(ordered=True)

    weights_outdir = os.path.join(outdir, "complete_weights")
    os.makedirs(weights_outdir, exist_ok=True)

    for view, df_W in W_complete_dict.items():
        # Transponemos (.T) para que: Filas = Genes/Variables, Columnas = Factores
        # Esto es lo ideal para GSEA y lectura manual
        df_W_final = df_W.T

        file_name = f"complete_weights_{view}_K{k}.csv"
        df_W_final.to_csv(os.path.join(weights_outdir, file_name))

    print(f"Archivos CSV de matrices completas guardados en: {outdir}")

    print(f"Finalizado K={k}\n")

print("Programa terminado")
