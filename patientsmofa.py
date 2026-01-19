import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Configuración de rutas
K_FOLDER = "K12"  # Cambia esto al modelo que hayas elegido
INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}"
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}"
GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"

# 1. Cargar datos
df_Z = pd.read_csv(os.path.join(INPUT_PATH, f"complete_factors_Z_{K_FOLDER}.csv"), index_col=0)
df_gv = pd.read_csv(GV_FILE, index_col=0)

# 2. Limpiar: solo factores reales para el heatmap
df_Z_only = df_Z[[c for c in df_Z.columns if "Factor" in c]]

# 3. Alinear pacientes
common_samples = df_Z_only.index.intersection(df_gv.index)
df_Z_plot = df_Z_only.loc[common_samples]
df_gv_plot = df_gv.loc[common_samples]

# 4. Crear el Clustermap
# Esto agrupará a los pacientes con perfiles biológicos similares
g = sns.clustermap(df_Z_plot.T, 
                   cmap="RdBu_r", 
                   center=0,
                   col_colors=df_gv_plot[['PREVIOUS_POLYPS', 'CRC_IN_FAMILY']], # Añade barras de color clínico
                   figsize=(12, 10),
                   xticklabels=False)

g.fig.suptitle(f'Z Matrix Heatmap - Patient Clusters ({K_FOLDER})', fontsize=16)
plt.savefig(os.path.join(OUTPUT_PATH, "z_matrix_heatmap.png"), dpi=300)
print(f"Heatmap guardado en {OUTPUT_PATH}")
