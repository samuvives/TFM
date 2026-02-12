import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. CONFIGURACIÓN DE RUTAS Y PARÁMETROS
# ==========================================
INPUT_PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
BASE_OUTPUT_PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/APPROACHSUBSETS/MOFAFLEX_FINAL_ANALYSIS"
# Ruta específica para las matrices Z (Factores)
Z_PATH_TEMPLATE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K{k}/complete_factors_Z_K{k}.csv"

K_LIST = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

# Mapeo: {Nombre_Vista_en_Output: Nombre_Archivo_en_Input}
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
# 2. FUNCIÓN DE CÁLCULO (RMSE DE RECONSTRUCCIÓN)
# ==========================================
def calcular_rmse_reconstruccion(path_input, path_weights, df_z):
    """
    Calcula el RMSE comparando el Input original con el producto Z * W^T
    """
    # Cargar datos (Input suele ser TSV, Weights suele ser CSV)
    df_in = pd.read_csv(path_input, sep='\t', index_col=0)
    df_w = pd.read_csv(path_weights, index_col=0)
    
    # 1. Alineación de factores entre Z y W (por si acaso difieren en orden o cantidad)
    common_factors = df_z.columns.intersection(df_w.columns)
    z_matrix = df_z[common_factors]
    w_matrix = df_w[common_factors]

    # 2. Reconstrucción matricial: Y_pred = Z x W.T
    y_pred = np.dot(z_matrix.values, w_matrix.values.T)
    df_pred = pd.DataFrame(y_pred, index=z_matrix.index, columns=df_w.index)

    # 3. Alineación de Muestras y Features
    # Solo comparamos lo que el modelo realmente entrenó (muestras y genes comunes)
    muestras_comunes = df_in.index.intersection(df_pred.index)
    features_comunes = df_in.columns.intersection(df_pred.columns)
    
    y_true_final = df_in.loc[muestras_comunes, features_comunes].values
    y_pred_final = df_pred.loc[muestras_comunes, features_comunes].values

    # 4. Cálculo matemático del RMSE
    # Usamos nanmean para manejar posibles valores faltantes en los datos crudos
    squared_error = (y_true_final - y_pred_final) ** 2
    rmse = np.sqrt(np.nanmean(squared_error))
    
    return rmse

# ==========================================
# 3. BUCLE DE PROCESAMIENTO
# ==========================================
results = []

print("Iniciando análisis de RMSE por factor...")

for k in K_LIST:
    print(f"--- Procesando K = {k} ---")
    weights_dir = f"{BASE_OUTPUT_PATH}/K{k}/complete_weights"
    path_z = Z_PATH_TEMPLATE.format(k=k)

    # Verificar existencia de archivos clave
    if not os.path.exists(weights_dir):
        print(f"  Aviso: No existe directorio de pesos para K{k}. Saltando...")
        continue
    if not os.path.exists(path_z):
        print(f"  Aviso: No existe archivo Z para K{k} en la ruta especificada. Saltando...")
        continue

    # Cargar Z una vez por cada K
    try:
        df_z = pd.read_csv(path_z, index_col=0)
    except Exception as e:
        print(f"  Error cargando Z en K{k}: {e}")
        continue

    # Procesar cada vista definida en el mapeo
    for view_name, input_filename in view_mapping.items():
        target_output = f"complete_weights_{view_name}_K{k}.csv"
        path_in = os.path.join(INPUT_PATH, input_filename)
        path_w = os.path.join(weights_dir, target_output)

        if os.path.exists(path_w) and os.path.exists(path_in):
            try:
                valor_rmse = calcular_rmse_reconstruccion(path_in, path_w, df_z)
                results.append({
                    'K': k,
                    'View': view_name,
                    'RMSE': valor_rmse
                })
            except Exception as e:
                print(f"  Error procesando vista {view_name} en K{k}: {e}")
        else:
            # Opcional: imprimir qué archivo falta
            pass

# Crear DataFrame con resultados
df_final = pd.DataFrame(results)

# ==========================================
# 4. GENERACIÓN Y GUARDADO DE RESULTADOS
# ==========================================
if not df_final.empty:
    # 1. Configuración de estilo del gráfico
    plt.figure(figsize=(14, 8))
    sns.set_style("whitegrid")
    sns.set_context("talk")

    # 2. Dibujar líneas de convergencia
    plot = sns.lineplot(
        data=df_final, 
        x='K', 
        y='RMSE', 
        hue='View', 
        marker='o', 
        linewidth=2.5,
        palette='tab20'
    )

    plt.title('Análisis de Plateau de RMSE por Vista Multi-ómica', pad=20)
    plt.xlabel('Número de Factores (K)')
    plt.ylabel('RMSE (Error de Reconstrucción)')
    plt.legend(title='Vistas', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()

    # 3. Guardar imagen en BASE_OUTPUT_PATH
    os.makedirs(BASE_OUTPUT_PATH, exist_ok=True)
    img_save_path = os.path.join(BASE_OUTPUT_PATH, 'analisis_plateau_rmse_mofa.png')
    plt.savefig(img_save_path, dpi=300)
    
    # 4. Guardar tabla CSV de respaldo
    csv_save_path = os.path.join(BASE_OUTPUT_PATH, 'metricas_rmse_resumen.csv')
    df_final.to_csv(csv_save_path, index=False)

    print(f"\nProceso finalizado con éxito.")
    print(f"Imagen guardada en: {img_save_path}")
    print(f"Datos guardados en: {csv_save_path}")

    # Mostrar el Delta RMSE (para identificar el codo numéricamente)
    print("\nVariación de RMSE (Delta) entre Ks sucesivos:")
    df_pivot = df_final.pivot(index='K', columns='View', values='RMSE')
    print(df_pivot.diff().round(4))

else:
    print("\nError: No se generaron datos. Revisa las rutas de los archivos Input y Output.")
