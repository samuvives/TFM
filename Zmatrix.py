# version of the Z heatmap with YES/NO in the guided factors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS
# ==========================================
K_FOLDER = "K12"
INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}"
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}"
GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"

os.makedirs(OUTPUT_PATH, exist_ok=True)

# ==========================================
# 2. CARGA Y ALINEACIÓN DE DATOS
# ==========================================
df_Z = pd.read_csv(os.path.join(INPUT_PATH, f"complete_factors_Z_{K_FOLDER}.csv"), index_col=0)
print("Z is:")
print(df_Z.head(3))
df_Z.index.name = 'sample_ID'
print("Z is:")
print(df_Z.head(3))

df_gv = pd.read_csv(GV_FILE, index_col=0)

# ==========================================
# 3. LABELING DE DATOS CLINICOS
# ==========================================
dfgvlabels = df_gv[['PREVIOUS_POLYPS', 'CRC_IN_FAMILY']].replace({1: "YES", 0: "NO"})
new_cols = [f"Factor {i}" for i in range(1, 13)]
dfgvlabels[new_cols] = ""
df_Z = df_Z.reindex(columns=dfgvlabels.columns) # is this necessary?

print("Z is:")
print(df_Z.head(3))
print("GV is:")
print(dfgvlabels.head(3))

# ==========================================
# 4. CREACIÓN DEL CLUSTERMAP CON NOMBRES
# ==========================================
# Ajustamos el tamaño de la figura para que quepan los nombres abajo
g = sns.clustermap(df_Z, 
                   cmap="RdBu_r", 
                   center=0,
                   annot=dfgvlabels,
                   fmt="",
                   figsize=(16, 12), # Un poco más ancho para los nombres
                   xticklabels=True,  # <--- AHORA ACTIVADO
                   yticklabels=True,
                   cbar_kws={'label': 'Factor Z-score'})

# Rotar los nombres de los pacientes para que sean legibles
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
# Ajustar los nombres de los factores
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

g.fig.suptitle(f'Z-Matrix Heatmap: Patient Stratification (K={K_FOLDER})', fontsize=18, y=1.05)

# ==========================================
# 6. SAVE RESULT
# ==========================================
save_path = os.path.join(OUTPUT_PATH, "z_matrix_heatmap_with_names.png")
# bbox_inches='tight' es vital para que no se corten los nombres de abajo al guardar
g.savefig(save_path, dpi=300, bbox_inches='tight')

print(f"Éxito: Heatmap con nombres guardado en {save_path}")
