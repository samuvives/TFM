import mofaflex as mfl
import pandas as pd
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS
# ==========================================
K_FOLDER = "K30"
BASE_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}"
MODEL_FILE = os.path.join(BASE_PATH, f"model_{K_FOLDER}.pkl")
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}"

os.makedirs(OUTPUT_PATH, exist_ok=True)

# ==========================================
# 2. CARGA DEL MODELO Y PROCESAMIENTO
# ==========================================
if os.path.exists(MODEL_FILE):
    print(f"Cargando modelo {K_FOLDER}...")
    model = mfl.MOFAFLEX.load(MODEL_FILE)
    
    # Obtener pesos (W)
    weights = model.get_weights(ordered=True)
    
    # Calculamos la actividad de cada vista por factor (Suma de cuadrados de pesos)
    # Esto indica cuánto contribuye cada ómica a cada factor
    view_activity_dict = {view: (W**2).sum(axis=1) for view, W in weights.items()}
    df_activity_all = pd.DataFrame(view_activity_dict)
    
    # Ordenar los factores por actividad total (de mayor a menor)
    total_activity = df_activity_all.sum(axis=1).sort_values(ascending=False)
    df_plot = df_activity_all.loc[total_activity.index]

    # ==========================================
    # 3. GENERACIÓN DEL GRÁFICO STACKED
    # ==========================================
    ax = df_plot.plot(kind='bar', 
                      stacked=True, 
                      figsize=(12, 7), 
                      colormap="tab20", 
                      edgecolor='white', 
                      linewidth=0.5)

    plt.title(f"Factor composition per view. K=30)", fontsize=16, pad=20)
    plt.ylabel("Sum of squared weights (Activity)", fontsize=12)
    plt.xlabel("Factores", fontsize=12)
    plt.xticks(rotation=45)

    plt.legend(title="Omic views)", 
               bbox_to_anchor=(1.05, 1), 
               loc='upper left', 
               fontsize=10, 
               frameon=False)

    plt.tight_layout()

    save_path = os.path.join(OUTPUT_PATH, "factor_composition_stacked_v2.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print(f"--- GRÁFICO GENERADO ---")
    print(f"Archivo guardado en: {save_path}")

else:
    print(f"Error: No se encontró el archivo {MODEL_FILE}")
