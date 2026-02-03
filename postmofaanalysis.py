import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# =======================================================
# 1. CONFIGURACIÓN DE RUTAS (BSC Project)
# =======================================================
INPUT_BASE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS"
OUTPUT_BASE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis"
GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"

os.makedirs(OUTPUT_BASE, exist_ok=True)

# Detectar modelos
k_folders = [d for d in os.listdir(INPUT_BASE) if d.startswith('K') and os.path.isdir(os.path.join(INPUT_BASE, d))]
df_gv = pd.read_csv(GV_FILE, index_col=0)
df_gv.index = df_gv.index.astype(str)
gv_names = df_gv.columns.tolist()

# =======================================================
# 2. BUCLE DE PROCESAMIENTO
# =======================================================
for k_dir in k_folders:
    print(f"\n>>> Analizando {k_dir}...")
    
    current_input_path = os.path.join(INPUT_BASE, k_dir)
    weights_dir = os.path.join(current_input_path, "complete_weights")
    factors_file = os.path.join(current_input_path, f"complete_factors_Z_{k_dir}.csv")
    
    current_output_path = os.path.join(OUTPUT_BASE, k_dir)
    os.makedirs(current_output_path, exist_ok=True)

    if not os.path.exists(factors_file):
        continue

    # --- CARGA Y LIMPIEZA DE COLUMNAS (Evita el ValueError) ---
    df_Z_raw = pd.read_csv(factors_file, index_col=0)
    df_Z_raw.index = df_Z_raw.index.astype(str)

    # Añadimos un prefijo "Model_" a todas las columnas que vienen de MOFA
    # Así PREVIOUS_POLYPS se convierte en Model_PREVIOUS_POLYPS y no choca con la clínica
    df_Z_raw.columns = [f"Model_{c}" for c in df_Z_raw.columns]

    common_samples = df_Z_raw.index.intersection(df_gv.index)
    df_combined = df_Z_raw.loc[common_samples].join(df_gv.loc[common_samples])

    # Definimos cuáles son los factores (las columnas que acabamos de renombrar)
    model_factors = df_Z_raw.columns.tolist()

    # --- 3.1.5: PATIENTS AND FACTORS ---
    # Usamos los dos primeros factores del modelo
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df_combined, x=model_factors[0], y=model_factors[1], 
                    hue=gv_names[0], palette="viridis", alpha=0.7)
    plt.title(f"Patient Map ({k_dir})\nX={model_factors[0]}, Y={model_factors[1]}")
    plt.savefig(os.path.join(current_output_path, "patient_map_Z1_Z2.png"))
    plt.close()

    # --- 3.1.7: RELATIONS WITH GUIDING VARIABLES ---
    # Calculamos correlación entre Factores del Modelo y Variables Clínicas Reales
    full_corr_mat = pd.DataFrame(index=model_factors, columns=gv_names)
    for gv in gv_names:
        for factor in model_factors:
            full_corr_mat.loc[factor, gv] = df_combined[factor].corr(df_combined[gv], method='spearman')

    for gv in gv_names:
        # Mejores 3 factores del modelo para esta variable clínica
        top_3 = full_corr_mat[gv].abs().sort_values(ascending=False).head(3).index
        for i, factor in enumerate(top_3):
            val = full_corr_mat.loc[factor, gv]
            plt.figure(figsize=(6, 6))
            sns.boxplot(data=df_combined, x=gv, y=factor, palette="Set2", showfliers=False)
            sns.stripplot(data=df_combined, x=gv, y=factor, color=".3", alpha=0.4)
            plt.title(f"{gv} vs {factor} (Rank {i+1})\nSpearman Rho: {val:.2f}")
            plt.tight_layout()
            plt.savefig(os.path.join(current_output_path, f"boxplot_{gv}_top{i+1}_{factor}.png"))
            plt.close()

    # --- 3.1.9: BIOLOGICAL RELEVANCE (Density) ---
    if os.path.exists(weights_dir):
        plt.figure(figsize=(10, 6))
        for f in os.listdir(weights_dir):
            if f.endswith(".csv"):
                view_name = f.replace(f"complete_weights_", "").replace(f"_{k_dir}.csv", "")
                W_df = pd.read_csv(os.path.join(weights_dir, f), index_col=0)
                sns.kdeplot(W_df.values.flatten(), label=view_name, bw_adjust=0.5)
        
        plt.title(f"Weight Density ({k_dir})")
        plt.xlim(-0.6, 0.6)
        plt.axvline(0, color='black', linestyle='--', alpha=0.3)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(current_output_path, "omics_density_dilution.png"))
        plt.close()

print(f"Post-análisis completado en: {OUTPUT_BASE}")
