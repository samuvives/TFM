import pandas as pd
import os

# 1. Define aquí tus archivos y sus sufijos correspondientes
# Estructura: "nombre_archivo.tsv": "SUFIJO"
configuracion = {
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_microbiota_SOFA_case_tumoral.tsv": "MIC",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_Metabolomics_data_case.tsv": "MET",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/renamed_lipidomics_data_case.tsv": "LIP"
}

for archivo, sufijo in configuracion.items():
    # Verificamos si el archivo existe en la carpeta para evitar errores
    if os.path.exists(archivo):
        try:
            # Leer el archivo TSV
            df = pd.read_csv(archivo, sep='\t')
            
            # 2. Renombrar columnas
            # La primera siempre será sample_ID, el resto llevará el sufijo del diccionario
            nuevas_cols = [
                "sample_ID" if i == 0 else f"{col}_{sufijo}" 
                for i, col in enumerate(df.columns)
            ]
            
            df.columns = nuevas_cols
            
            # 3. Guardar el archivo
            df.to_csv(archivo, sep='\t', index=False)
            print(f"Procesado: {archivo} (Sufijo: _{sufijo})")
            
        except Exception as e:
            print(f"Error al procesar {archivo}: {e}")
    else:
        print(f"El archivo '{archivo}' no se encontró en el directorio.")
