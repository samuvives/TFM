import pandas as pd
import os
import sys

# Arrange all the files in the same order

# 1. Files
archivos = [
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/VC/MPA_GT_1_2.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/EXPRESSION/tpmexpression.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/renamed_Metabolomics_data_case.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/renamed_lipidomics_data_case.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/renamed_microbiota_SOFA_case_tumoral.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DEL_Patient_Gene_Matrix.tsv", 
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_DUP_Patient_Gene_Matrix.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INS_Patient_Gene_Matrix.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_INV_Patient_Gene_Matrix.tsv",
    "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV/SV_TRA_Patient_Gene_Matrix.tsv"
] 

columna_id = "sample_ID"

print("Starting alignment", flush=True)

# Stablish reference order
try:
    print(f"Reference IDs from the file: {os.path.basename(archivos[0])}", flush=True)
    orden_maestro = pd.read_csv(archivos[0], sep='\t', usecols=[columna_id])[columna_id].tolist()
    print(f"Order stablished, {len(orden_maestro)} samples.", flush=True)
except Exception as e:
    print(f"CRITICAL ERROR reading reference file: {e}", flush=True)
    sys.exit(1)

# Process and overwrite the files
for i, ruta_archivo in enumerate(archivos, 1):
    nombre_archivo = os.path.basename(ruta_archivo)
    print(f"[{i}/{len(archivos)}] Processing: {nombre_archivo}...", end=" ", flush=True)
    
    if os.path.exists(ruta_archivo):
        try:
            # MEMORY OPTIMIZATION: 
            # engine='c' es más rápido para archivos grandes.
            # low_memory=False evita avisos de tipos de datos mixtos en matrices grandes.
            df = pd.read_csv(ruta_archivo, sep='\t', low_memory=False, engine='c')
            
            # Rearrange
            # put the sample_ID as index, apply the order, and get sample_ID back as a column
            df_reordenado = df.set_index(columna_id).reindex(orden_maestro).reset_index()
            
            # Overwrite the original file
            # Na_rep='' asegura que los valores vacíos se guarden como celdas vacías (estándar TSV)
            df_reordenado.to_csv(ruta_archivo, sep='\t', index=False, na_rep='NA')
            
            print("DONE", flush=True)
            del df
            del df_reordenado
            
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
    else:
        print(f"ERROR: File not found in the path.", flush=True)

print("\n Process finished with success.", flush=True)
