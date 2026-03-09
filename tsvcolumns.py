# 
import pandas as pd
import os

config = {
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DEL_Patient_Gene_Matrix.tsv": "DEL",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DUP_Patient_Gene_Matrix.tsv": "DUP",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INS_Patient_Gene_Matrix.tsv": "INS",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INV_Patient_Gene_Matrix.tsv": "INV",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_TRA_Patient_Gene_Matrix.tsv": "TRA"
}

# new
from tsvcolumns2 import addsufixes_tsv
for filepath, suffix in config.items():
    addsufixes_tsv(filepath, suffix)

# old
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
