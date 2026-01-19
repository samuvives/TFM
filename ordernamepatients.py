import pandas as pd
import os
import sys

# 1. Configuración 
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

print(">>> Iniciando proceso de alineación optimizado...", flush=True)

# 2. Establecer el ORDEN DE REFERENCIA de forma ligera
# Optimizamos leyendo SOLO la columna sample_ID para no cargar 400MB innecesariamente
try:
    print(f">>> Leyendo IDs de referencia desde: {os.path.basename(archivos[0])}", flush=True)
    orden_maestro = pd.read_csv(archivos[0], sep='\t', usecols=[columna_id])[columna_id].tolist()
    print(f">>> Orden maestro establecido con {len(orden_maestro)} muestras.", flush=True)
except Exception as e:
    print(f"ERROR CRÍTICO al leer el archivo de referencia: {e}", flush=True)
    sys.exit(1)

# 3. Procesar y sobreescribir
for i, ruta_archivo in enumerate(archivos, 1):
    nombre_archivo = os.path.basename(ruta_archivo)
    print(f"[{i}/{len(archivos)}] Procesando: {nombre_archivo}...", end=" ", flush=True)
    
    if os.path.exists(ruta_archivo):
        try:
            # OPTIMIZACIÓN DE MEMORIA: 
            # engine='c' es más rápido para archivos grandes.
            # low_memory=False evita avisos de tipos de datos mixtos en matrices grandes.
            df = pd.read_csv(ruta_archivo, sep='\t', low_memory=False, engine='c')
            
            # Reordenar
            # reindex es ideal para MOFA porque si falta un ID, crea la fila con NaNs
            df_reordenado = df.set_index(columna_id).reindex(orden_maestro).reset_index()
            
            # Sobreescribir el archivo original
            # Na_rep='' asegura que los valores vacíos se guarden como celdas vacías (estándar TSV)
            df_reordenado.to_csv(ruta_archivo, sep='\t', index=False, na_rep='NA')
            
            print("HECHO ✅", flush=True)
            
            # Liberar memoria explícitamente (buena práctica en HPC)
            del df
            del df_reordenado
            
        except Exception as e:
            print(f"ERROR: {e}", flush=True)
    else:
        print(f"ERROR: Archivo no encontrado en la ruta.", flush=True)

print("\n>>> Proceso finalizado con éxito.", flush=True)
