import pandas as pd
import os

# --- RUTAS ---
folder_path = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/GV'
excel_path = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/CORRESPONDENCIA.xlsx'
reference_tsv = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/SV_TRA_Patient_Gene_Matrix.tsv'
ids_to_remove = [70350514, 5024299, 166088, 5506597]

# 1. Cargar el Excel y forzar ambas columnas a string
df_map = pd.read_excel(excel_path)
df_map['GV'] = df_map['GV'].astype(str)
df_map['NEWGV'] = df_map['NEWGV'].astype(str)

# Crear el diccionario de mapeo
name_mapping = dict(zip(df_map['GV'], df_map['NEWGV']))

# 2. Obtener el orden deseado desde el TSV (convertido a string)
df_ref = pd.read_csv(reference_tsv, sep='\t')
desired_order = df_ref['sample_ID'].astype(str).tolist()

# 3. Procesar archivos CSV
for filename in os.listdir(folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(folder_path, filename)
        
        # Cargar CSV
        df = pd.read_csv(file_path)
        
        # Renombrar columna
        df = df.rename(columns={'sampleID': 'sample_ID'})
        
        # Convertir IDs actuales a string para la comparación y eliminación
        df['sample_ID'] = df['sample_ID'].astype(str)
        
        # Eliminar filas prohibidas (usando strings)
        ids_to_remove_str = [str(x) for x in ids_to_remove]
        df = df[~df['sample_ID'].isin(ids_to_remove_str)]
        
        # Aplicar el cambio a nombres nuevos (que ya son strings en el mapping)
        df['sample_ID'] = df['sample_ID'].replace(name_mapping)
        
        # 4. Reordenar según el TSV de referencia
        # Ponemos la columna como índice para usar el reindexado rápido
        df = df.set_index('sample_ID')
        
        # Reindexamos: esto asegura el orden del TSV y descarta lo que no esté en él
        # Usamos intersection para evitar errores si el TSV pide un ID que no existe en el CSV
        existing_order = [x for x in desired_order if x in df.index]
        df = df.reindex(existing_order).reset_index()
        
        # Guardar cambios
        df.to_csv(file_path, index=False)
        print(f"Procesado correctamente: {filename}")

print("--- Proceso finalizado con éxito ---")
