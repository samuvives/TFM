import pandas as pd
import os
# pasamos de csv a tsv y los nombres de columnas dejan de tener espacios para tener guiones
# Definimos la ruta absoluta
directorio = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC'

# Obtener la lista de archivos
archivos = os.listdir(directorio)

for nombre_archivo in archivos:
    if nombre_archivo.endswith('.csv'):
        # Construimos la ruta completa de entrada
        ruta_entrada = os.path.join(directorio, nombre_archivo)
        print(f"Procesando: {ruta_entrada}...")
        
        try:
            # 1. Leer el CSV usando la ruta completa
            df = pd.read_csv(ruta_entrada)
            
            # 2. Limpiar nombres de columnas
            df.columns = [col.strip().replace(' ', '_').replace('"', '') for col in df.columns]
            
            # 3. Crear el nombre de salida (.tsv) en la misma ruta
            nombre_tsv = nombre_archivo.rsplit('.', 1)[0] + '.tsv'
            ruta_salida = os.path.join(directorio, nombre_tsv)
            
            # 4. Guardar como TSV
            df.to_csv(ruta_salida, sep='\t', index=False)
            
            print(f"Guardado en: {ruta_salida}")
            
        except Exception as e:
            print(f"Error en {nombre_archivo}: {e}")

print("\nProceso finalizado con éxit")
