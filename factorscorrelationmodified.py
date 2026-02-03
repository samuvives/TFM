# study the correlation between factors
# guided and non-guided factors
# heatmap with guided and non-guided
# study the non-guided to help stablish the adequate number of factors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os


def corrfactorshist(df, namecolumn, hist_title, save_path):
    # 1. Definir el ancho que deseas (ejemplo: 0.1)
    ancho_barra = 0.05
    limites = np.arange(-1, 1 + ancho_barra, ancho_barra)
    plt.figure(figsize=(10, 6))
    sns.histplot(df[namecolumn], bins=limites, kde=True, color='skyblue', edgecolor='black')

    plt.title('Correlation values histogram (-1 to 1)')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.xticks(limites, rotation=45, ha="right", rotation_mode="anchor")
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def processcorrmatrix(corr_matrix):
    # 1. Limpiar espacios en los nombres para evitar errores visuales
    corr_matrix.columns = corr_matrix.columns.str.replace(' ', '')
    corr_matrix.index = corr_matrix.index.str.replace(' ', '')
    
    # 2. Obtener solo el triángulo superior sin la diagonal (k=1)
    # Esto elimina las correlaciones repetidas (A-B y B-A) y las de un factor consigo mismo
    mask = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
    
    # 3. "Apilar" la matriz para pasar de tabla a lista de pares
    # stack() elimina automáticamente los valores nulos (el triángulo inferior que ocultamos)
    bardf = corr_matrix.where(mask).stack().reset_index()
    bardf.columns = ['F1', 'F2', 'corrbarvalues']
    
    # 4. Crear el nombre del par para el eje X
    bardf['corrbarnames'] = bardf['F1'] + " vs " + bardf['F2']

    return bardf

def corrfactorsbarplot(bardf, hmtitle, save_path):
    """Does a barplot of the correlation values"""

    # --- Diseño del Gráfico ---
    plt.figure(figsize=(14, 7))
    colores = ['#4C72B0' if v > 0 else '#C44E52' for v in bardf["corrbarvalues"]]
    
    sns.barplot(x="corrbarnames", y="corrbarvalues", data=bardf, palette=colores)
    
    plt.axhline(0, color='black', linewidth=1)
    plt.ylim(-1.1, 1.1)
    plt.xticks(rotation=90)
    plt.xlabel("Pares de Factores")
    plt.ylabel("Coeficiente de Spearman")
    plt.title(hmtitle, fontsize=16, pad=20)
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"--- BARPLOT GENERADO ---")
    print(f"Archivo guardado en: {save_path}")


def corrfactorsheatmap(corr_matrix, hmtitle, save_path):
    """Takes a matrix and creates the heatmap"""
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix, 
                annot=True,
                fmt=".2f",
                cmap="RdBu_r",
                center=0,
                vmin=-1, vmax=1,
                square=True,
                linewidths=.5,
                cbar_kws={"shrink": .8, "label": "Spearman Correlation"})

    # Formatear etiquetas para que no se corten
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.title(hmtitle, fontsize=16, pad=20)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"--- HEATMAP GENERADO ---")
    print(f"Archivo guardado en: {save_path}")


K_NUMBER = "K12"
INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_NUMBER}"
OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_NUMBER}"
os.makedirs(OUTPUT_PATH, exist_ok=True)

# Load factors file
file_z = os.path.join(INPUT_PATH, f"complete_factors_Z_K12.csv")
df_Z = pd.read_csv(file_z, index_col=0)

# Calculamos la correlación de todos los factores (guiados y latentes)
corr_matrix = df_Z.corr(method='spearman')

# extraemos aparte solo los no guiados
corr_matrix_ng = corr_matrix.filter(regex="^Factor")
corr_matrix_ng = corr_matrix_ng.filter(regex="^Factor", axis=0)

# creamos heatmap con guiados y no guiados
guided_nonguided_hmtitle = f'Inter-Factor Correlation Heatmap ({K_NUMBER})\nGuided & Non-Guided Factors'
save_path_heatmap = os.path.join(OUTPUT_PATH, "factor_correlation_heatmap.png")
corrfactorsheatmap(corr_matrix, guided_nonguided_hmtitle, save_path_heatmap)

# procesamos para mas visualizaciones
valuesdf = processcorrmatrix(corr_matrix_ng)

# creamos barplot
nonguided_bartitle = f"Correlation values ({K_NUMBER})\nNon-Guided Factors"
save_path_barplot = os.path.join(OUTPUT_PATH, "factor_correlation_barplot.png")
corrfactorsbarplot(valuesdf, nonguided_bartitle, save_path_barplot)

# creamos histograma
nonguided_histtitle = 'Correlation values histogram (-1 to 1)'
save_path_hist = os.path.join(OUTPUT_PATH, "factor_correlation_hist.png")
corrfactorshist(valuesdf, "corrbarvalues", nonguided_histtitle, save_path_hist)
