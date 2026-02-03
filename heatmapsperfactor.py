import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS
# ==========================================
K_FOLDER = "K12"
INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}"
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}/individual_factors_large_text"

os.makedirs(OUTPUT_PATH, exist_ok=True)

# ==========================================
# 2. CARGA DE DATOS
# ==========================================
file_z = os.path.join(INPUT_PATH, f"complete_factors_Z_{K_FOLDER}.csv")
df_Z = pd.read_csv(file_z, index_col=0)

# ==========================================
# 3. GENERACIÓN DE HEATMAPS (TEXTO GRANDE)
# ==========================================
for factor in df_Z.columns:
    data_factor = df_Z[[factor]].T
    
    # Aumentamos el ancho (25) y el alto (6) para dar aire a las etiquetas grandes
    plt.figure(figsize=(25, 6))
    
    g = sns.clustermap(data_factor,
                       cmap="RdBu_r",
                       center=0,
                       col_cluster=True,
                       row_cluster=False,
                       xticklabels=True, 
                       yticklabels=[factor],
                       cbar_pos=(0.005, 0.4, 0.01, 0.3),
                       figsize=(25, 7))

    # AJUSTE DE TEXTO: Tamaño 9 y negrita suave para mejor lectura
    plt.setp(g.ax_heatmap.get_xticklabels(), 
             rotation=90, 
             fontsize=12, 
             weight='normal')
    
    # Ajuste del nombre del Factor (Eje Y)
    plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=14, weight='bold')
    
    # Limpiar etiquetas de ejes
    g.ax_heatmap.set_xlabel("Patient Identifiers", fontsize=12, labelpad=15)
    g.ax_heatmap.set_ylabel("")
    
    # Título estilizado
    g.fig.suptitle(f'Factor Analysis: {factor} - Patient Clustering', 
                   fontsize=20, 
                   y=1.05)
    
    # Guardar con margen extra para nombres largos
    clean_name = factor.replace(".", "_")
    save_file = os.path.join(OUTPUT_PATH, f"Large_Detail_{clean_name}.png")
    
    # dpi=300 garantiza que al ampliar el PDF/Imagen el nombre no se pixele
    g.savefig(save_file, dpi=300, bbox_inches='tight')
    
    plt.close()

print(f"--- PROCESO COMPLETADO ---")
print(f"Heatmaps con nombres grandes guardados en: {OUTPUT_PATH}")
