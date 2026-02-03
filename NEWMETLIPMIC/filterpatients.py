import pandas as pd
import os

# Ruta del directorio
directorio = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC'

# Lista de IDs a eliminar (como strings para evitar problemas de formato)
muestras_a_eliminar = ['70350514', '5024299', '166088', '5506597']

# Obtener la lista de archivos TSV
archivos = [f for f in os.listdir(directorio) if f.endswith('.csv')]

for nombre_archivo in archivos:
    ruta_completa = os.path.join(directorio, nombre_archivo)
    print(f"Filtrando muestras en: {nombre_archivo}...")
    
    try:
        # 1. Leer el TSV
        df = pd.read_csv(ruta_completa)
        
        if 'Sample' in df.columns:
            # 2. Asegurar que la columna Sample se trata como string para comparar correctamente
            filas_antes = len(df)
            
            # 3. Filtrar: mantener solo las filas cuyo 'Sample' NO esté en nuestra lista
            df['Sample'] = df['Sample'].astype(str)
            df = df[~df['Sample'].isin(muestras_a_eliminar)]
            
            filas_despues = len(df)
            eliminadas = filas_antes - filas_despues
            
            # 4. Sobrescribir el archivo con los datos filtrados
            df.to_csv(ruta_completa, sep='\t', index=False)
            print(f"Hecho. Se eliminaron {eliminadas} filas.")
        else:
            print(f"Omitido: La columna 'Sample' no existe en {nombre_archivo}")
            
    except Exception as e:
        print(f"Error procesando {nombre_archivo}: {e}")

print("\nFiltrado completado.")
