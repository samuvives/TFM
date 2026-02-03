# guided and non-guided factors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS
# ==========================================
K_FOLDER = "K12"
INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}"
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}"

os.makedirs(OUTPUT_PATH, exist_ok=True)

# ==========================================
# 2. CARGA DE DATOS
# ==========================================
file_z = os.path.join(INPUT_PATH, f"complete_factors_Z_K12.csv")
df_Z = pd.read_csv(file_z, index_col=0)

# ==========================================
# 3. CÁLCULO DE CORRELACIÓN (Spearman)
# ==========================================
# Calculamos la correlación de todos los factores (guiados y latentes)
corr_matrix = df_Z.corr(method='spearman')

# ==========================================
# 4. GENERACIÓN DEL HEATMAP (Sin Clustering)
# ==========================================
plt.figure(figsize=(12, 10))

# Usamos sns.heatmap directamente en lugar de clustermap
sns.heatmap(corr_matrix, 
            annot=True,             # Mostrar los valores
            fmt=".2f",              # Dos decimales
            cmap="RdBu_r",          # Rojo (pos), Azul (neg)
            center=0,               # Blanco en 0
            vmin=-1, vmax=1,        # Rango completo de correlación
            square=True,            # Celdas cuadradas
            linewidths=.5,          # Separación entre celdas
            cbar_kws={"shrink": .8, "label": "Spearman Correlation"})

# Formatear etiquetas para que no se corten
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.title(f'Inter-Factor Correlation Heatmap (K={K_FOLDER})\nGuided & Latent Factors', fontsize=16, pad=20)

# ==========================================
# 5. GUARDAR Y RESUMIR
# ==========================================
save_path = os.path.join(OUTPUT_PATH, "factor_correlation_flat.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')

print(f"--- HEATMAP GENERADO ---")
print(f"Archivo guardado en: {save_path}")

# Identificar las correlaciones más significativas entre factores guiados y latentes
# Filtramos la matriz para ver solo filas de guiados vs columnas de latentes
guided_vars = ['PREVIOUS_POLYPS', 'CRC_IN_FAMILY'] # Ajusta si tienes más
latent_factors = [c for c in df_Z.columns if "Factor" in c]

if all(v in corr_matrix.index for v in guided_vars):
    print("\nAsociaciones clave (Guiados vs Latentes):")
    summary = corr_matrix.loc[guided_vars, latent_factors]
    # Mostrar las 3 correlaciones más altas por cada variable guiada
    for var in guided_vars:
        top_3 = summary.loc[var].abs().sort_values(ascending=False).head(3)
        print(f"\nTop factores asociados a {var}:")
        for f, val in top_3.items():
            real_val = summary.loc[var, f]
            print(f"  -> {f}: {real_val:.3f}")
