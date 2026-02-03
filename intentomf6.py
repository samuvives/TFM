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
import dask.dataframe as dd # Importamos dask para asegurar compatibilidad

warnings.filterwarnings("ignore", category=FutureWarning)

# parameters
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproach"
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
N_FACTORS_LIST = [20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]

# Configuración de vista
VIEWS_CONFIG = {
    "Metabolomics": {"path": f'{INPUT_DIR}/renamed_Metabolomics_data_case.tsv', "likelihood": "Normal", "scale": True},
    "Lipidomics": {"path": f'{INPUT_DIR}/renamed_lipidomics_data_case.tsv', "likelihood": "Normal", "scale": True},
    "Microbiota": {"path": f'{INPUT_DIR}/renamed_microbiota_SOFA_case_tumoral.tsv', "likelihood": "Normal", "scale": True},
    "SV_DEL": {"path": f'{INPUT_DIR}/SV_DEL_Patient_Gene_Matrix.tsv', "likelihood": "Bernoulli", "scale": False},
    "SV_DUP": {"path": f'{INPUT_DIR}/SV_DUP_Patient_Gene_Matrix.tsv', "likelihood": "Bernoulli", "scale": False},
    "SV_INS": {"path": f'{INPUT_DIR}/SV_INS_Patient_Gene_Matrix.tsv', "likelihood": "Bernoulli", "scale": False},
    "SV_INV": {"path": f'{INPUT_DIR}/SV_INV_Patient_Gene_Matrix.tsv', "likelihood": "Bernoulli", "scale": False},
    "SV_TRA": {"path": f'{INPUT_DIR}/SV_TRA_Patient_Gene_Matrix.tsv', "likelihood": "Bernoulli", "scale": False},
    "VC_01": {"path": f'{INPUT_DIR}/MPA_GT_0_1.tsv', "likelihood": "Bernoulli", "scale": False},
    "VC_11": {"path": f'{INPUT_DIR}/MPA_GT_1_1.tsv', "likelihood": "Bernoulli", "scale": False},
    "VC_12": {"path": f'{INPUT_DIR}/MPA_GT_1_2.tsv', "likelihood": "Bernoulli", "scale": False},
    "EXPRESSION": {"path": f'{INPUT_DIR}/tpmexpression.tsv', "likelihood": "Normal", "scale": True}
}

GV_FILE = f"{INPUT_DIR}/matrizGV4_mapeado.csv"
base_output = os.path.join(WORKING_DIR, "Single_Scenario_Analysis_GV1")
os.makedirs(base_output, exist_ok=True)

# =======================================================
# 2. FUNCIONES DE DIAGNÓSTICO
# =======================================================

def run_diagnostics_single(model, gv_names, run_tag, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    try:
        mfl.pl.variance_explained(model)
        plt.title(f'Varianza - {run_tag}')
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_variance.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    try:
        factors = model.get_factors()
        gv_cols = [c for c in factors.columns if c in gv_names]
        lf_cols = [c for c in factors.columns if c not in gv_names]
        factors_reordered = factors[gv_cols + lf_cols]
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(factors_reordered, cmap="RdBu_r", center=0)
        if gv_cols:
            plt.axvline(x=len(gv_cols), color='black', linewidth=2)
        plt.title(f'Pacientes vs Factores - {run_tag}')
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_factors_Z.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    for view in model.data.views:
        try:
            mfl.pl.top_loadings(model, views=view, n_features=20)
            plt.savefig(os.path.join(fig_dir, f"{run_tag}_loadings_{view}.png"), bbox_inches='tight')
            plt.close()
        except: plt.close()

# =======================================================
# 3. EJECUCION
# =======================================================

# 3.1 Carga de datos optimizada
print(">>> Cargando modalidades y aplicando escalado si es necesario...")
adata_dict = {}
for mod, config in VIEWS_CONFIG.items():
    print(f"    -> Procesando: {mod}")
    df = pd.read_csv(config["path"], sep="\t", index_col=0).sort_index()
    
    if config["scale"]:
        scaler = StandardScaler()
        df = pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns)
    
    # Creamos AnnData asegurando float32 para ahorrar memoria
    adata_dict[mod] = anndata.AnnData(
        X=df.values.astype(np.float32), 
        obs=pd.DataFrame(index=df.index), 
        var=pd.DataFrame(index=df.columns)
    )

mdata = MuData(adata_dict)

# 3.2 Carga de GV y preparación de MuData (UNA SOLA VEZ)
print(">>> Inyectando variables guía (GV)...")
gv_df = pd.read_csv(GV_FILE, index_col=0).sort_index()
gv_names = gv_df.columns.tolist()

for mod in mdata.mod.keys():
    mdata.mod[mod].obs = mdata.mod[mod].obs.join(gv_df)

# 3.3 Bucle de entrenamiento por K
for k in N_FACTORS_LIST:
    run_tag = f"AllModalities_GV1_K{k}"
    run_dir = os.path.join(base_output, f"K{k}")
    os.makedirs(run_dir, exist_ok=True)
    
    print(f"\n--- Iniciando Entrenamiento: K={k} ---")
    
    try:
        # Configuración del modelo (usamos mdata directamente al haber inyectado GV antes)
        model = mfl.MOFAFLEX(
            mdata,
            mfl.ModelOptions(
                likelihoods={m: VIEWS_CONFIG[m]["likelihood"] for m in VIEWS_CONFIG.keys()},
                guiding_vars_likelihoods={c: "Bernoulli" for c in gv_names},
                n_factors=k
            ),
            mfl.TrainingOptions(seed=42),
            mfl.DataOptions(guiding_vars_obs_keys=gv_names)
        )
        
        print(f">>> Generando diagnósticos para K={k}...")
        run_diagnostics_single(model, gv_names, run_tag, run_dir)
        
        print(f">>> Guardando pesos para K={k}...")
        weights = model.get_weights()
        for view, W in weights.items():
            W.to_csv(os.path.join(run_dir, f"weights_{view}_K{k}.csv"))
            
        print(f"✅ Completado con éxito: K={k}")
        
    except Exception as e:
        print(f"❌ Error en la ejecución de K={k}: {e}")

print(f"\n🚀 Proceso finalizado. Resultados en: {base_output}")
