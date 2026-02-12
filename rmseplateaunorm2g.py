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

# Diccionario de likelihoods (Bernoulli para vistas SV/VC, Gaussian para el resto)
view_likelihood = {}
for view in view_mapping.keys():
    if view.startswith("SV_") or view.startswith("VC_"):
        view_likelihood[view] = "bernoulli"
    else:
        view_likelihood[view] = "gaussian"

print("Distribuciones por vista:")
for view, likelihood in view_likelihood.items():
    print(f"  {view}: {likelihood}")

# ==========================================
# 2. FUNCIÓN DE CÁLCULO (nRMSE) con soporte para Bernoulli
# ==========================================
def sigmoid(x):
    x = np.clip(x, -50, 50)  # numerical stability
    return 1.0 / (1.0 + np.exp(-x))

def calcular_metricas_reconstruccion(path_input, path_weights, df_z, likelihood="gaussian"):
    # Load data
    df_in = pd.read_csv(path_input, sep="\t", index_col=0)
    df_w  = pd.read_csv(path_weights, index_col=0)

    # Align factors
    common_factors = df_z.columns.intersection(df_w.columns)
    z_matrix = df_z[common_factors]
    w_matrix = df_w[common_factors]

    # Reconstruction on linear predictor
    lin_pred = np.dot(z_matrix.values, w_matrix.values.T)  # (samples, features)

    # Apply link if needed
    if likelihood.lower() == "bernoulli":
        y_pred = sigmoid(lin_pred)     # probabilities in [0,1]
    else:
        y_pred = lin_pred              # gaussian mean

    df_pred = pd.DataFrame(y_pred, index=z_matrix.index, columns=df_w.index)

    # Align samples/features
    m_comunes = df_in.index.intersection(df_pred.index)
    f_comunes = df_in.columns.intersection(df_pred.columns)

    y_true = df_in.loc[m_comunes, f_comunes].values.astype(float)
    y_p    = df_pred.loc[m_comunes, f_comunes].values.astype(float)

    # RMSE (NaN-safe)
    diff = y_true - y_p
    rmse = float(np.sqrt(np.nanmean(diff * diff)))

    # Normalization
    if likelihood.lower() == "bernoulli":
        # Null predictor = prevalence per feature (or global). Here: per-feature.
        p_null = np.nanmean(y_true, axis=0, keepdims=True)  # (1, features)
        rmse_null = float(np.sqrt(np.nanmean((y_true - p_null) ** 2)))
        nrmse = rmse / rmse_null if rmse_null > 0 else rmse
    else:
        std_dev = float(np.nanstd(y_true))
        nrmse = rmse / std_dev if std_dev > 0 else rmse

    return rmse, nrmse

# ==========================================
# 3. BUCLE PRINCIPAL
# ==========================================
results = []

for k in K_LIST:
    print(f"\nCalculando para K={k}...")
    path_z = Z_PATH_TEMPLATE.format(k=k)
    weights_dir = f"{BASE_OUTPUT_PATH}/K{k}/complete_weights"

    if not os.path.exists(path_z) or not os.path.exists(weights_dir):
        print(f"  Advertencia: No se encontraron datos para K={k}")
        continue

    df_z = pd.read_csv(path_z, index_col=0)
    views_procesadas = 0

    for view, inp_file in view_mapping.items():
        path_in = os.path.join(INPUT_PATH, inp_file)
        path_w = os.path.join(weights_dir, f"complete_weights_{view}_K{k}.csv")

        if os.path.exists(path_w):
            try:
                # Obtener el tipo de distribución para esta vista
                likelihood = view_likelihood[view]
                rmse_val, nrmse_val = calcular_metricas_reconstruccion(
                    path_in, path_w, df_z, likelihood=likelihood
                )
                results.append({
                    'K': k, 
                    'View': view, 
                    'Likelihood': likelihood,
                    'RMSE_Crudo': rmse_val, 
                    'nRMSE': nrmse_val
                })
                views_procesadas += 1
                print(f"  ✓ {view} ({likelihood}): nRMSE = {nrmse_val:.4f}")
            except Exception as e:
                print(f"  ✗ Error en {view} K{k}: {e}")
        else:
            print(f"  ✗ No se encontró archivo de pesos para {view} K{k}")

    print(f"  Procesadas: {views_procesadas} vistas para K={k}")

# Crear DataFrame con resultados
df_final = pd.DataFrame(results)

# ==========================================
# 4. GRÁFICO Y GUARDADO
# ==========================================
if not df_final.empty:
    # Configurar estilo
    plt.figure(figsize=(14, 8))
    sns.set_style("whitegrid")
    
    # Crear gráfico separando por tipo de distribución
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Subplot 1: Vistas Gaussianas
    gaussian_views = df_final[df_final['Likelihood'] == 'gaussian']
    if not gaussian_views.empty:
        sns.lineplot(data=gaussian_views, x='K', y='nRMSE', hue='View', 
                    marker='o', linewidth=2, ax=axes[0])
        axes[0].set_title('Vistas Gaussianas (Continuas)', fontsize=14)
        axes[0].set_ylabel('nRMSE (Normalizado por Desviación Estándar)')
        axes[0].set_xlabel('Número de Factores (K)')
        axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Vista')
    
    # Subplot 2: Vistas Bernoulli
    bernoulli_views = df_final[df_final['Likelihood'] == 'bernoulli']
    if not bernoulli_views.empty:
        sns.lineplot(data=bernoulli_views, x='K', y='nRMSE', hue='View', 
                    marker='s', linewidth=2, ax=axes[1])
        axes[1].set_title('Vistas Bernoulli (SV/VC)', fontsize=14)
        axes[1].set_ylabel('nRMSE (Normalizado por Prevalencia)')
        axes[1].set_xlabel('Número de Factores (K)')
        axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Vista')
    
    plt.suptitle('Convergencia del Modelo por Tipo de Distribución', fontsize=16)
    plt.tight_layout()

    # Guardar en la ruta de salida
    save_img = os.path.join(BASE_OUTPUT_PATH, "plot_nrmse_convergencia.png")
    save_csv = os.path.join(BASE_OUTPUT_PATH, "metricas_completas_rmse.csv")
    
    plt.savefig(save_img, dpi=300, bbox_inches='tight')
    df_final.to_csv(save_csv, index=False)
    
    print(f"\n{'='*60}")
    print(f"Finalizado. Archivos guardados en:")
    print(f"  Imagen: {save_img}")
    print(f"  CSV: {save_csv}")
    print(f"{'='*60}")
    
    # Mostrar tabla resumen del nRMSE
    print("\nTabla Resumen nRMSE (Valores bajos = mejor ajuste):")
    summary_table = df_final.pivot_table(
        index=['View', 'Likelihood'], 
        columns='K', 
        values='nRMSE',
        aggfunc='first'
    ).round(4)
    print(summary_table)
    
    # Guardar también la tabla resumen
    summary_csv = os.path.join(BASE_OUTPUT_PATH, "resumen_nrmse_pivot.csv")
    summary_table.to_csv(summary_csv)
    print(f"\nTabla resumen guardada en: {summary_csv}")
    
else:
    print("\nNo se encontraron resultados. Verifique las rutas de entrada.")
