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
    df = corr_matrix.where(mask).stack().reset_index()
    df.columns = ['F1', 'F2', 'corrbarvalues']
    
    # 4. Crear el nombre del par para el eje X
    df['corrbarnames'] = df['F1'] + " vs " + df['F2']

    return df

klist = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
correlationsumdict = {}
correlationmeandict = {}
for k in klist:
    K_NUMBER = "K" + str(k)
    INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_NUMBER}"
    OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/histcorrelationanalysis"
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    # Load factors file
    file_z = os.path.join(INPUT_PATH, f"complete_factors_Z_{K_NUMBER}.csv")
    df_Z = pd.read_csv(file_z, index_col=0)

    # Calculamos la correlación de todos los factores (guiados y latentes)
    corr_matrix = df_Z.corr(method='spearman')

    # extraemos aparte solo los no guiados
    corr_matrix_ng = corr_matrix.filter(regex="^Factor")
    corr_matrix_ng = corr_matrix_ng.filter(regex="^Factor", axis=0)

    # procesamos para mas visualizaciones
    valuesdf = processcorrmatrix(corr_matrix_ng)

    # creamos histograma
    nonguided_histtitle = f"Correlation values histogram (-1 to 1) {K_NUMBER}"
    save_path_hist = os.path.join(OUTPUT_PATH, "factor_correlation_hist_{K_NUMBER}.png")
    corrfactorshist(valuesdf, "corrbarvalues", nonguided_histtitle, save_path_hist)

    # nos quedamos con la suma de los valores de correlacion
    total_suma = (valuesdf['corrbarvalues'] ** 2).sum()
    correlationsumdict[K_NUMBER] = total_suma

    # nos quedamos con la media de los valores de correlacion
    promedio = (valuesdf['corrbarvalues'] ** 2).mean()
    correlationmeandict[K_NUMBER] = promedio

# sum 
nombres = list(correlationsumdict.keys())
valores = list(correlationsumdict.values())

plt.scatter(nombres, valores)

plt.title('Sum value of factor correlation per number of factors')
plt.xlabel('Number of factors')
plt.ylabel('Sum of correlation values')
plt.savefig(os.path.join(OUTPUTPATH, "sumcorrelationineachnumberoffactors.png"))
plt.close()


# mean
nombres = list(correlationmeandict.keys())
valores = list(correlationmeandict.values())

plt.scatter(nombres, valores)

plt.title('Mean value of factor correlation per number of factors')
plt.xlabel('Number of factors')
plt.ylabel('Mean correlation value')
plt.savefig(os.path.join(OUTPUTPATH, "meancorrelationineachnumberoffactors.png"))
plt.close()
