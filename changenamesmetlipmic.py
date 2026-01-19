import pandas as pd
import os
# cambiamos los nombres de pacientes a la referencia y quitamos los 4 cuya muestra de transcriptomica es defectuosa

# 1. CONFIGURACIÓN DE RUTAS
base_path = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC/'
path_mapeo = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/omicasmiquelprocessing/CORRESPONDENCIA.xlsx'

# Lista de archivos a procesar
archivos_a_procesar = [
    'Metabolomics_data_case.tsv', 
    'lipidomics_data_case.tsv', 
    'microbiota_SOFA_case_tumoral.tsv'
]

# 2. CARGAR MAPEO Y CREAR DICCIONARIO
df_mapeo = pd.read_excel(path_mapeo)
dict_renombrar = dict(zip(
    df_mapeo['MET'].astype(str).str.strip(), 
    df_mapeo['NEWMET'].astype(str).str.strip()
))

patientstoeliminate = ["70350514", "5024299", "166088", "5506597"]

# 3. BUCLE DE PROCESAMIENTO
for nombre_archivo in archivos_a_procesar:
    path_entrada = os.path.join(base_path, nombre_archivo)
    
    # Verificar si el archivo existe antes de abrirlo
    if os.path.exists(path_entrada):
        print(f"Procesando: {nombre_archivo}...")
        
        # Cargar archivo (asumiendo que el ID es la primera columna y la usamos como índice)
        df = pd.read_csv(path_entrada, sep="\t", index_col=0)

        # Renombrar las filas (índice)
        df.rename(index=dict_renombrar, inplace=True)

        # Eliminar pacientes específicos del índice
        df.drop(index=patientstoeliminate, inplace=True, errors='ignore')

        # Guardar el resultado (sobreescribiendo o creando un prefijo 'proc_')
        path_salida = os.path.join(base_path, f"renamed_{nombre_archivo}")
        df.to_csv(path_salida)
        
        print(f"  - Guardado en: {path_salida}")
        print(f"  - Dimensiones finales: {df.shape}")
    else:
        print(f"Advertencia: El archivo {nombre_archivo} no se encontró en la ruta.")

print("\nProceso completado para todos los archivos.")
