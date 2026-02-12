import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS Y PARÁMETROS
# ==========================================
INPUT_PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
BASE_OUTPUT_PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS"
Z_PATH_TEMPLATE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K{k}/complete_factors_Z_K{k}.csv"

K_LIST = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

view_mapping = {
    "EXPRESSION": "tpmexpression.tsv",
    "Microbiota": "renamed_microbiota_SOFA_case_tumoral.tsv",
    "Metabolomics": "renamed_Metabolomics_data_case.tsv",
    "Lipidomics": "renamed_lipidomics_data_case.tsv",
    "SV_DEL": "SV_DEL_Patient_Gene_Matrix.tsv",
    "SV_DUP": "SV_DUP_Patient_Gene_Matrix.tsv",
    "SV_INS": "SV_INS_Patient_Gene_Matrix.tsv",
    "SV_INV": "SV_INV_Patient_Gene_Matrix.tsv",
    "SV_TRA": "SV_TRA_Patient_Gene_Matrix.tsv",
    "VC_11": "MPA_GT_1_1.tsv",
    "VC_12": "MPA_GT_1_2.tsv",
    "VC_01": "MPA_GT_0_1.tsv"
}

# ==========================================
# 2. FUNCIÓN DE CÁLCULO (nRMSE)
# ==========================================
def calcular_nrmse_reconstruccion(path_input, path_weights, df_z):
    # Cargar datos
    df_in = pd.read_csv(path_input, sep='\t', index_col=0)
    df_w = pd.read_csv(path_weights, index_col=0)
    
    # Alineación de factores
    common_factors = df_z.columns.intersection(df_w.columns)
    z_matrix = df_z[common_factors]
    w_matrix = df_w[common_factors]

    # Reconstrucción: Y_pred = Z * W.T
    y_pred = np.dot(z_matrix.values, w_matrix.values.T)
    df_pred = pd.DataFrame(y_pred, index=z_matrix.index, columns=df_w.index)

    # Alineación de muestras y variables
    m_comunes = df_in.index.intersection(df_pred.index)
    f_comunes = df_in.columns.intersection(df_pred.columns)
    
    y_true = df_in.loc[m_comunes, f_comunes].values
    y_p = df_pred.loc[m_comunes, f_comunes].values

    # RMSE crudo
    rmse = np.sqrt(np.nanmean((y_true - y_p)**2))
    
    # NORMALIZACIÓN: nRMSE = RMSE / Desviación Estándar de los datos originales
    # Esto permite que todas las vistas compitan en la misma escala (0 a 1 aprox)
    std_dev = np.nanstd(y_true)
    nrmse = rmse / std_dev if std_dev != 0 else rmse
    
    return rmse, nrmse

# ==========================================
# 3. BUCLE PRINCIPAL
# ==========================================
results = []

for k in K_LIST:
    print(f"Calculando para K={k}...")
    path_z = Z_PATH_TEMPLATE.format(k=k)
    weights_dir = f"{BASE_OUTPUT_PATH}/K{k}/complete_weights"

    if not os.path.exists(path_z) or not os.path.exists(weights_dir):
        continue

    df_z = pd.read_csv(path_z, index_col=0)

    for view, inp_file in view_mapping.items():
        path_in = os.path.join(INPUT_PATH, inp_file)
        path_w = os.path.join(weights_dir, f"complete_weights_{view}_K{k}.csv")

        if os.path.exists(path_w):
            try:
                rmse_val, nrmse_val = calcular_nrmse_reconstruccion(path_in, path_w, df_z)
                results.append({
                    'K': k, 
                    'View': view, 
                    'RMSE_Crudo': rmse_val, 
                    'nRMSE': nrmse_val
                })
            except Exception as e:
                print(f"Error en {view} K{k}: {e}")

df_final = pd.DataFrame(results)

# ==========================================
# 4. GRÁFICO Y GUARDADO
# ==========================================
if not df_final.empty:
    plt.figure(figsize=(12, 7))
    sns.set_style("whitegrid")
    
    # Graficamos el nRMSE (Normalizado) para que todas las líneas sean visibles
    sns.lineplot(data=df_final, x='K', y='nRMSE', hue='View', marker='o', linewidth=2)
    
    plt.title('Convergencia del Modelo (nRMSE por Vista)', fontsize=15)
    plt.ylabel('nRMSE (Normalizado por Desviación Estándar)')
    plt.xlabel('Número de Factores (K)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Guardar en la ruta de salida
    save_img = os.path.join(BASE_OUTPUT_PATH, "plot_nrmse_convergencia.png")
    save_csv = os.path.join(BASE_OUTPUT_PATH, "metricas_completas_rmse.csv")
    
    plt.savefig(save_img, dpi=300)
    df_final.to_csv(save_csv, index=False)
    
    print(f"\nFinalizado. Archivos guardados en {BASE_OUTPUT_PATH}")
    
    # Mostrar tabla resumen del nRMSE
    print("\nTabla Resumen nRMSE (Valores bajos = mejor ajuste):")
    print(df_final.pivot(index='View', columns='K', values='nRMSE').round(4))
