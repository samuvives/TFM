import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- CONFIGURACIÓN DE RUTAS ---
BASE_DIR = "/home/vant/Escritorio/TFM/resultadosmf/Single_Scenario_Analysis_GV1"
K_FOLDERS = [f"K{i}" for i in range(20, 42, 2)] # Genera K20, K22... K40

print(f">>> Iniciando generación de gráficas en: {BASE_DIR}")

for k_name in K_FOLDERS:
    k_path = os.path.join(BASE_DIR, k_name)
    
    if not os.path.exists(k_path):
        print(f"--- [!] Saltando {k_name}: No se encuentra la carpeta ---")
        continue

    print(f"--- Procesando {k_name} ---")
    
    # Crear carpeta para las nuevas gráficas manuales
    output_figs = os.path.join(k_path, "manual_loadings")
    os.makedirs(output_figs, exist_ok=True)

    # Buscar archivos de pesos
    weight_files = [f for f in os.listdir(k_path) if f.startswith("weights_") and f.endswith(".csv")]

    for file in weight_files:
        # Extraer nombre de la vista (ej. Metabolomics)
        view_name = file.replace("weights_", "").split("_K")[0]
        
        try:
            # Leer el CSV
            df = pd.read_csv(os.path.join(k_path, file), index_col=0)
            
            # Si los factores están en el índice, transponemos para tener Genes/Metabolitos en filas
            if "Factor 1" in df.index:
                df = df.T

            # Graficamos el Factor 1 (suele ser el más importante)
            factor = "Factor 1"
            if factor in df.columns:
                plt.figure(figsize=(8, 10))
                
                # Seleccionar los 30 mayores por valor absoluto
                top_30 = df[factor].abs().sort_values(ascending=False).head(30).index
                data_plot = df.loc[top_30, factor].sort_values()
                
                # Colores: Azul positivo, Rojo negativo
                colors = ['steelblue' if x > 0 else 'indianred' for x in data_plot]
                
                data_plot.plot(kind='barh', color=colors)
                plt.title(f"Top 30 Loadings: {view_name}\n({k_name} - {factor})", fontsize=12)
                plt.xlabel("Weight Value")
                plt.axvline(0, color='black', lw=1, ls='-')
                plt.grid(axis='x', linestyle='--', alpha=0.5)
                
                # Guardar
                plt.savefig(os.path.join(output_figs, f"{view_name}_{factor}.png"), bbox_inches='tight')
                plt.close()
                
        except Exception as e:
            print(f"    [!] Error en {file}: {e}")

print(f"\n🚀 ¡Listo! Revisa las carpetas 'manual_loadings' dentro de cada K.")
