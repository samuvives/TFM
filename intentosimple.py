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

# parameters
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproach"
INPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/hugoapproach"
N_FACTORS_LIST = [20, 22, 24, 26, 28, 30]

VIEWS_CONFIG = {
    "Metabolomics": {"path": f'{INPUT_DIR}/metfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Met"},
    "Lipidomics": {"path": f'{INPUT_DIR}/lipfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Lip"},
    "Microbiota": {"path": f'{INPUT_DIR}/micfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Mic"},
    "SV": {"path": f'{INPUT_DIR}/SVfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "SV"},
    "AEB": {"path": f'{INPUT_DIR}/AEfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "AEB"},
    "AS": {"path": f'{INPUT_DIR}/ASfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "AS"}
}

GV_FILE = f"{INPUT_DIR}/GV1filtered.csv"
base_output = os.path.join(WORKING_DIR, "Single_Scenario_Analysis_GV1")
os.makedirs(base_output, exist_ok=True)

# =======================================================
# 2. FUNCIONES DE DIAGNÓSTICO
# =======================================================

def run_diagnostics_single(model, gv_names, run_tag, output_dir):
    """Genera las gráficas solicitadas para cada run de K."""
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    # 1. Varianza Explicada
    try:
        mfl.pl.variance_explained(model)
        plt.title(f'Varianza - {run_tag}')
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_variance.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    # 2. Heatmap de Factores vs Pacientes (Matriz Z) - REORDENADO
    try:
        factors = model.get_factors()
        gv_cols = [c for c in factors.columns if c in gv_names]
        lf_cols = [c for c in factors.columns if c not in gv_names]
        factors_reordered = factors[gv_cols + lf_cols]
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(factors_reordered, cmap="RdBu_r", center=0)
        if gv_cols:
            plt.axvline(x=len(gv_cols), color='black', linewidth=2)
        plt.title(f'Pacientes vs Factores (GVs a la izquierda) - {run_tag}')
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_factors_Z.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    # 3. Top Loadings por Vista
    for view in model.data.views:
        try:
            mfl.pl.top_loadings(model, views=view, n_features=20)
            plt.savefig(os.path.join(fig_dir, f"{run_tag}_loadings_{view}.png"), bbox_inches='tight')
            plt.close()
        except: plt.close()

# =======================================================
# 3. EJECUCION (Bucle sobre K)
# =======================================================

# 3.1 Carga de datos inicial
print(">>> Cargando modalidades con prefijos...")
current_modalities = {}
for mod, config in VIEWS_CONFIG.items():
    df = pd.read_csv(config["path"], index_col=0).sort_index()
    df.columns = [f"{config['prefix']}_{c}" for c in df.columns]
    if config["scale"]:
        df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)
    current_modalities[mod] = df

adata_dict = {name: anndata.AnnData(X=df.values.astype(float), 
                                   obs=pd.DataFrame(index=df.index), 
                                   var=pd.DataFrame(index=df.columns)) 
              for name, df in current_modalities.items()}
mdata = MuData(adata_dict)

# 3.2 Carga de GV de family history
gv_df = pd.read_csv(GV_FILE, index_col=0).sort_index()
gv_names = gv_df.columns.tolist()

# 3.3 Bucle de entrenamiento por K
for k in N_FACTORS_LIST:
    run_tag = f"AllModalities_GV1_K{k}"
    run_dir = os.path.join(base_output, f"K{k}")
    os.makedirs(run_dir, exist_ok=True)
    
    print(f"\n--- Iniciando Run: K={k} ---")
    
    # Inyectar GV en el objeto MuData
    mdata_run = mdata.copy()
    for mod in VIEWS_CONFIG.keys():
        mdata_run.mod[mod].obs = mdata_run.mod[mod].obs.join(gv_df)
    
    try:
        # Configuración del modelo
        model = mfl.MOFAFLEX(
            mdata_run,
            mfl.ModelOptions(
                likelihoods={m: VIEWS_CONFIG[m]["likelihood"] for m in VIEWS_CONFIG.keys()},
                guiding_vars_likelihoods={c: "Bernoulli" for c in gv_names},
                n_factors=k
            ),
            mfl.TrainingOptions(seed=42),
            mfl.DataOptions(guiding_vars_obs_keys=gv_names)
        )
        
        # Generar todas las gráficas
        run_diagnostics_single(model, gv_names, run_tag, run_dir)
        
        # Guardar pesos
        weights = model.get_weights()
        for view, W in weights.items():
            W.to_csv(os.path.join(run_dir, f"weights_{view}_K{k}.csv"))
            
        print(f"Completado con éxito: K={k}")
        
    except Exception as e:
        print(f"Error en la ejecución de K={k}: {e}")

print(f"\nProceso finalizado. Todas las gráficas están en: {base_output}")
