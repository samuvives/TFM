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
from matplotlib.patches import Patch

warnings.filterwarnings("ignore", category=FutureWarning)

# =======================================================
# 1. CONFIGURACIÓN Y RUTAS
# =======================================================
WORKING_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/totalapproach"
INPUT_DIR = "/home/bsc/bsc236340/MOFAINPUT/hugoapproach"
N_FACTORS_LIST = [20, 22, 24, 26, 28, 30]

VIEWS_CONFIG = {
    "Metabolomics": {"path": f'{INPUT_DIR}/metfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Met"},
    "Lipidomics": {"path": f'{INPUT_DIR}/lipfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Lip"},
    "Microbiota": {"path": f'{INPUT_DIR}/micfiltered.csv', "likelihood": "Normal", "scale": True, "prefix": "Mic"},
    "SV": {"path": f'{INPUT_DIR}/SVfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "SV"},
    "AEB": {"path": f'{INPUT_DIR}/AEfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "AEB"},
    "AS": {"path": f'{INPUT_DIR}/ASfiltered.csv', "likelihood": "Bernoulli", "scale": False, "prefix": "AS"}
}

guiding_files = [f"{INPUT_DIR}/GV{i}filtered.csv" for i in range(1, 7)]
base_output = os.path.join(WORKING_DIR, "Results_Total_Analysis")
os.makedirs(base_output, exist_ok=True)

# =======================================================
# 2. FUNCIONES DE APOYO
# =======================================================

def clean_feature_name(name):
    for cfg in VIEWS_CONFIG.values():
        prefix = f"{cfg['prefix']}_"
        if name.startswith(prefix):
            return name.replace(prefix, "", 1)
    return name

def run_full_diagnostics(model, gv_names, run_tag, output_dir):
    fig_dir = os.path.join(output_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    # 1. Varianza Explicada
    try:
        mfl.pl.variance_explained(model)
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_variance.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    # 2. Heatmap de Factores vs Pacientes (Matriz Z)
    try:
        factors = model.get_factors()
        gv_cols = [c for c in factors.columns if c in gv_names]
        lf_cols = [c for c in factors.columns if c not in gv_names]
        factors_reordered = factors[gv_cols + lf_cols]
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(factors_reordered, cmap="RdBu_r", center=0)
        if gv_cols:
            plt.axvline(x=len(gv_cols), color='black', linewidth=2)
        plt.title(f'Z Matrix (GVs Left) - {run_tag}')
        plt.savefig(os.path.join(fig_dir, f"{run_tag}_factors_Z.png"), bbox_inches='tight')
        plt.close()
    except: plt.close()

    # 3. Top Loadings por Vista (Acceso corregido a model.views)
    # Usamos model.views en lugar de model.data.views
    views_to_plot = getattr(model, 'views', [])
    for view in views_to_plot:
        try:
            mfl.pl.top_loadings(model, views=view, n_features=15)
            plt.savefig(os.path.join(fig_dir, f"{run_tag}_loadings_{view}.png"), bbox_inches='tight')
            plt.close()
        except: plt.close()

# =======================================================
# 3. PROCESAMIENTO
# =======================================================

full_traceability_data = []
unique_features_pool = set()

scenarios = {"All_Modalities": list(VIEWS_CONFIG.keys())}
for mod in VIEWS_CONFIG.keys():
    scenarios[f"Excluding_{mod}"] = [m for m in VIEWS_CONFIG.keys() if m != mod]

for scenario_name, active_views in scenarios.items():
    print(f"\n>>> ESCENARIO: {scenario_name}")
    scenario_accumulator = []
    
    current_modalities = {}
    for mod in active_views:
        df = pd.read_csv(VIEWS_CONFIG[mod]["path"], index_col=0).sort_index()
        df.columns = [f"{VIEWS_CONFIG[mod]['prefix']}_{c}" for c in df.columns]
        if VIEWS_CONFIG[mod]["scale"]:
            df = pd.DataFrame(StandardScaler().fit_transform(df), index=df.index, columns=df.columns)
        current_modalities[mod] = df

    adata_dict = {name: anndata.AnnData(X=df.values.astype(float), 
                                       obs=pd.DataFrame(index=df.index), 
                                       var=pd.DataFrame(index=df.columns)) 
                  for name, df in current_modalities.items()}
    mdata = MuData(adata_dict)

    for i, gv_path in enumerate(guiding_files, 1):
        gv_df = pd.read_csv(gv_path, index_col=0).sort_index()
        gv_names = gv_df.columns.tolist()
        gv_id = f"GV{i}"
        
        for k in N_FACTORS_LIST:
            run_tag = f"{scenario_name}_{gv_id}_K{k}"
            run_dir = os.path.join(base_output, scenario_name, gv_id, f"K{k}")
            os.makedirs(run_dir, exist_ok=True)
            print(f"  -> {run_tag}")

            mdata_run = mdata.copy()
            for mod in active_views: 
                mdata_run.mod[mod].obs = mdata_run.mod[mod].obs.join(gv_df)
            
            try:
                model = mfl.MOFAFLEX(
                    mdata_run, 
                    mfl.ModelOptions(likelihoods={m: VIEWS_CONFIG[m]["likelihood"] for m in active_views},
                                     guiding_vars_likelihoods={c: "Bernoulli" for c in gv_names},
                                     n_factors=k),
                    mfl.TrainingOptions(seed=42),
                    mfl.DataOptions(guiding_vars_obs_keys=gv_names)
                )
                
                run_full_diagnostics(model, gv_names, run_tag, run_dir)
                
                weights = model.get_weights()
                for view, W in weights.items():
                    if W is not None:
                        active_cols = [c for c in W.columns if c in gv_names]
                        for f_name in active_cols:
                            df_l = pd.DataFrame({
                                "Feature_Prefixed": W.index,
                                "Feature_Base": [clean_feature_name(n) for n in W.index],
                                "Loading_Abs": W[f_name].abs(),
                                "Factor_Name": f_name,
                                "View": view,
                                "GV_Set": gv_id,
                                "Scenario": scenario_name,
                                "K": k
                            })
                            scenario_accumulator.append(df_l)
                            full_traceability_data.append(df_l)
                            
            except Exception as e:
                print(f"      Error: {e}")

    if scenario_accumulator:
        scen_df = pd.concat(scenario_accumulator)
        top_50 = scen_df.groupby("Feature_Base")["Loading_Abs"].max().sort_values(ascending=False).head(50).index.tolist()
        unique_features_pool.update(top_50)

# =======================================================
# 4. OUTPUTS FINALES
# =======================================================
if full_traceability_data:
    pd.concat(full_traceability_data).to_excel(os.path.join(base_output, "FINAL_TRACEABILITY.xlsx"), index=False)

if unique_features_pool:
    all_dfs = []
    for cfg in VIEWS_CONFIG.values():
        df = pd.read_csv(cfg["path"], index_col=0).sort_index()
        all_dfs.append(df)
    full_df = pd.concat(all_dfs, axis=1)
    
    # Media de columnas con el mismo nombre (sin prefijo)
    full_df_clean = full_df.groupby(by=full_df.columns, axis=1).mean()
    
    heatmap_feats = [f for f in unique_features_pool if f in full_df_clean.columns]
    h_plot = full_df_clean[heatmap_feats].T
    h_scaled = (h_plot.sub(h_plot.mean(axis=1), axis=0)).div(h_plot.std(axis=1), axis=0)

    g = sns.clustermap(h_scaled, cmap="RdBu_r", center=0, figsize=(15, 12), yticklabels=True, xticklabels=False)
    plt.savefig(os.path.join(base_output, "GLOBAL_POOL_HEATMAP.png"), dpi=300, bbox_inches="tight")

print(f"\nProceso finalizado. Resultados en: {base_output}")
