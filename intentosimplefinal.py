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

def zfactors(file, outdir, gv_names):
    """Aplica en seccion 3.1, Function to know the inactive factors"""
    fdata = pd.read_csv(file, index_col=0)
    fdata = fdata[fdata.index.str.contains("Factor")]
    
    max_valor = fdata.iloc[:, 0].max()

    # threshold
    thr = 0.01 * max_valor

    mask = fdata.iloc[:, 0] < thr
    inactive_list = fdata.index[mask].tolist()
    if len(inactive_list) == 0:
        print("All factors exceed the activity threshold")
    else:
        print(f"The inactive factors are : {', '.join(inactive_list)}")


# =======================================================
# 1. PARAMETERS
# =======================================================
print("[1] Initializing parameters")
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
# 2. FILE UPLOAD 
# =======================================================
print("[2] Loading views")
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
# 3 MAIN LOOP
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

    # 3.1 Checking of Z activity
    Z_dict = model.get_factors(ordered=True)
    Z = Z_dict[list(Z_dict.keys())[0]]
    factor_activity = (Z ** 2).mean(axis=0)
    zfactorsfile = os.path.join(outdir, "factor_activity.csv")
    factor_activity.to_csv(zfactorsfile)
    zfactors(zfactorsfile, outdir, gv_names)
    print("3.1 done")

    # 3.2 W Weights & Stacked Bar of the distribution of Activity
    weights = model.get_weights(ordered=True) 
    view_activity_dict = {view: (W**2).sum(axis=1) for view, W in weights.items()}
    df_activity_all = pd.DataFrame(view_activity_dict)
    total_activity = df_activity_all.sum(axis=1).sort_values(ascending=False)
    
    df_activity_all.loc[total_activity.index].plot(kind='bar', stacked=True, figsize=(10, 5), colormap="tab10") # tab10 is too low
    plt.title(f"Composición de Factores (K={k})")
    plt.savefig(os.path.join(figdir, "factor_composition_stacked.png"))
    plt.close()
    print("3.2 done")

    # ===================================================
    # 3.10 GLOBAL RANKING
    # ===================================================
    print(f"[3.10] Calculating global ranking global for K={k}")
    all_feats = []
    for view, W in weights.items():
        imp = (W**2).sum(axis=0)
        df_imp = imp.to_frame(name='GlobalScore')
        df_imp['View'] = view
        df_imp['PctContrib'] = (df_imp['GlobalScore'] / df_imp['GlobalScore'].sum()) * 100
        all_feats.append(df_imp)

    # Definimos df_global uniendo todos los datos
    df_global = pd.concat(all_feats).sort_values(by='GlobalScore', ascending=False)
    print("todo okey en 3.10")

    # ===================================================
    # 3.11 LIDERES POR VISTA (TOP 3)
    # ===================================================
    print(f"[3.11] Extrayendo los Top 3 líderes de cada vista")

    # Reseteamos el índice para que el nombre del gen/metabolito sea una columna llamada 'Feature'
    df_global_reset = df_global.reset_index().rename(columns={'index': 'Feature'})

    # Extraemos el top 3 de cada grupo 'View'
    top_leaders = (
        df_global_reset.groupby('View')
        .apply(lambda x: x.sort_values('GlobalScore', ascending=False).head(3))
        .reset_index(drop=True)
    )

    top_leaders.to_csv(os.path.join(outdir, "top_3_leaders_per_view.csv"), index=False)
    print("todo okey en 3.11")

    # ===================================================
    # 3.12 HEATMAP DE LÍDERES (INTEGRADO)
    # ===================================================
    print("[3.12] Generando Heatmap de líderes")
    list_w = []
    for _, row in top_leaders.iterrows():
        v = row['View']
        f = row['Feature']

        # Accedemos a los pesos originales usando la vista y el nombre de la variable
        ws = weights[v][f]
        ws.name = f"{v}: {f}"
        list_w.append(ws)

    plt.figure(figsize=(12, 10))
    sns.heatmap(pd.concat(list_w, axis=1).T, cmap="RdBu_r", center=0)
    plt.title(f"Integrated Leaders Heatmap (K={k})")
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "integrated_leaders_heatmap.png"))
    plt.close()
    print("todo okey en 3.12")

    # ===================================================
    # 3.13 LOADING PLOTS (ENFOQUE DE LÍDERES) - BLINDADO
    # ===================================================
    print(f"[3.13] Generando Loading Plots para los líderes de cada vista (K={k})")

    # Calculamos el número de vistas para el layout
    n_views = len(weights)
    fig, axes = plt.subplots(n_views, 1, figsize=(10, 4 * n_views), sharex=False)

    # Ajuste por si solo hay una vista (axes no sería una lista)
    if n_views == 1:
        axes = [axes]

    for i, (view, W) in enumerate(weights.items()):
        # Buscamos el líder #1 de esta vista en nuestro DataFrame top_leaders
        try:
            # Filtramos el top_leaders para esta vista y tomamos el primero
            leader_row = top_leaders[top_leaders['View'] == view].iloc[0]
            leader_name = leader_row['Feature']

            # Extraemos los pesos del leader weights (Series: Factores -> Pesos)
            leader_loadings = W[leader_name]

            # Graficamos
            sns.barplot(x=leader_loadings.index, y=leader_loadings.values, ax=axes[i], palette="vlag")

            axes[i].set_title(f"Líder de {view}: {leader_name}", fontsize=12, fontweight='bold')
            axes[i].set_ylabel("Weight (W)")
            axes[i].axhline(0, color='black', linewidth=0.8)
            axes[i].tick_params(axis='x', rotation=45)
            print("todo okey en 3.13")

        except Exception as e:
            print(f"   [!] Error procesando líder para la vista {view}: {e}")
            axes[i].set_title(f"Error en vista {view}")

    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "leaders_loading_focus.png"), dpi=300)
    plt.close()

    # ===================================================
    # 3.14 CLINICAL CORRELATION (FACTORS VS GV)
    # ===================================================
    print("[3.14] Calculating clinical correlation")
    common_samples = Z.index.intersection(gv_df.index)
    clinical_corr = pd.DataFrame(index=Z.columns, columns=gv_names)
    for f_col in Z.columns:
        for g_col in gv_names:
            clinical_corr.loc[f_col, g_col] = Z.loc[common_samples, f_col].corr(gv_df.loc[common_samples, g_col], method='spearman')
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(clinical_corr.astype(float), cmap="PiYG", center=0, annot=True, fmt=".2f")
    plt.title(f"Clinical Correlation Heatmap (K={k})")
    plt.tight_layout()
    plt.savefig(os.path.join(figdir, "clinical_correlation_heatmap.png"))
    plt.close()
    print("3.14 done")


    # ===================================================
    # 3.15 EXPORT OF FULL MATRICES (Z AND W)
    # ===================================================
    print(f"[3.15] Saving factors and weights matrices (K={k})")

    # 1. FACTORS MATRIX (Z)
    Z_complete_dict = model.get_factors(ordered=True)
    group_name = list(Z_complete_dict.keys())[0]
    df_Z_final = Z_complete_dict[group_name]

    df_Z_final.to_csv(os.path.join(outdir, f"complete_factors_Z_K{k}.csv"))

    # 2. WEIGHTS MATRIX (W) PER VIEW
    W_complete_dict = model.get_weights(ordered=True)

    weights_outdir = os.path.join(outdir, "complete_weights")
    os.makedirs(weights_outdir, exist_ok=True)

    for view, df_W in W_complete_dict.items():
        # Transponemos (.T) para que: Filas = Genes/Variables, Columnas = Factores
        # Esto es lo ideal para GSEA y lectura manual
        df_W_final = df_W.T

        file_name = f"complete_weights_{view}_K{k}.csv"
        df_W_final.to_csv(os.path.join(weights_outdir, file_name))

    print(f"CSV files of full matrices saved in: {outdir}")
    print("3.15 done")

    print(f"Completed K={k}\n")

print("Program completed")
