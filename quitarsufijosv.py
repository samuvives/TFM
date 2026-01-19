import pandas as pd
import os

# Tu ruta exacta según el comando pwd
ruta = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/"
columna_id = "sample_ID"

# Cambiamos al directorio de trabajo para que pd.read_csv encuentre los archivos
os.chdir(ruta)

for archivo in os.listdir(ruta):
    if archivo.endswith(".tsv"):
        # 1. Leer archivo
        df = pd.read_csv(archivo, sep='\t')
        
        # 2. Modificar la columna eliminando el sufijo
        # Usamos .str.replace o .str.removesuffix
        df[columna_id] = df[columna_id].astype(str).str.replace("_sv", "", regex=False)
        
        # 3. Guardar (sobreescribiendo el original)
        # index=False es vital para no añadir una columna extra de números
        df.to_csv(archivo, sep='\t', index=False)
        
        print(f"Modificado: {archivo}")
