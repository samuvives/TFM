import pandas as pd
import os

# --- CONFIGURACIÓN ---
BASE_DIR = "/home/vant/Escritorio/TFM/resultadosmf/Single_Scenario_Analysis_GV1"
# Buscamos todas las carpetas que empiecen por K
K_FOLDERS = [d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d)) and d.startswith('K')]
TOP_N = 20  # Número de features top a extraer por cada modalidad/factor

print(f">>> Iniciando análisis multimodal en: {BASE_DIR}")

for k_name in sorted(K_FOLDERS, key=lambda x: int(x[1:])):
    k_path = os.path.join(BASE_DIR, k_name)
    output_file = os.path.join(k_path, f"Ranking_Global_Features_{k_name}.csv")
    
    # 1. Listar archivos de pesos en la carpeta actual
    weight_files = [f for f in os.listdir(k_path) if f.startswith("weights_") and f.endswith(".csv")]
    
    if not weight_files:
        print(f"--- [!] Saltando {k_name}: No se encontraron archivos de pesos ---")
        continue

    print(f"--- Procesando {k_name} ({len(weight_files)} modalidades) ---")
    resumen_k = []

    for file in weight_files:
        # Extraer nombre de la modalidad
        view_name = file.replace("weights_", "").split("_K")[0]
        
        try:
            df = pd.read_csv(os.path.join(k_path, file), index_col=0)
            
            # Asegurar orientación: Features en filas, Factores en columnas
            if "Factor 1" in df.index:
                df = df.T
            
            for factor in df.columns:
                # Extraer las top N features por valor absoluto para este factor y esta vista
                top_data = df[factor].abs().sort_values(ascending=False).head(TOP_N)
                
                for feature, abs_weight in top_data.items():
                    real_weight = df.loc[feature, factor]
                    resumen_k.append({
                        "Factor": factor,
                        "Modalidad": view_name,
                        "Feature": feature,
                        "Weight": real_weight,
                        "Abs_Weight": abs_weight
                    })
        except Exception as e:
            print(f"    [!] Error procesando {file} en {k_name}: {e}")

    # 2. Crear DataFrame de la carpeta K actual
    if resumen_k:
        df_k = pd.DataFrame(resumen_k)
        
        # Ordenar: primero por Factor, luego por importancia absoluta (independiente de la modalidad)
        df_k = df_k.sort_values(by=["Factor", "Abs_Weight"], ascending=[True, False])
        
        # 3. Guardar CSV
        df_k.to_csv(output_file, index=False)
        print(f"    ✅ Guardado ranking en: {output_file}")
    else:
        print(f"    [!] No se generaron datos para {k_name}")

print(f"\n🚀 Proceso finalizado para todas las carpetas.")
