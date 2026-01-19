import pandas as pd

# 1. CARGAR LOS ARCHIVOS
path_mapeo = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/CORRESPONDENCIA.xlsx'
path_tpm = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/EXPRESSIONPROCESSING/tpmssplicing/transcript_tpm_matrix_suppa.tsv'
path_salida = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/EXPRESSIONPROCESSING/tpmexpression.csv'

df_mapeo = pd.read_excel(path_mapeo)
df_principal = pd.read_csv(path_tpm, sep="\t")

# 2. CREAR EL DICCIONARIO DE CORRESPONDENCIA (Evitando errores de tipo string)
dict_renombrar = dict(zip(
    df_mapeo['EXPRESSION'].astype(str).str.strip(), 
    df_mapeo['NEWEXPRESSION'].astype(str).str.strip()
))

# 3. RENOMBRAR LAS COLUMNAS (Pacientes)
df_principal.rename(columns=dict_renombrar, inplace=True)

# 4. ELIMINAR PACIENTES ESPECIFICOS
patientstoeliminate = ["70350514", "5024299", "166088", "5506597"]
# errors='ignore' evita que el script falle si un ID no existe en las columnas
df_principal.drop(columns=patientstoeliminate, inplace=True, errors='ignore')

# 5. PREPARAR PARA TRANSPONER
# Identificamos la primera columna 'transcript_id'
nombre_columna_genes = df_principal.columns[0]

# Fijamos los genes como índice para que al transponer se conviertan en cabeceras
df_principal.set_index(nombre_columna_genes, inplace=True)

# 6. TRANSPONER
# Ahora: Filas = Pacientes, Columnas = Genes
df_final = df_principal.T

# 7. RECUPERAR EL ID DEL PACIENTE
# Al transponer, el ID del paciente quedó en el índice. Lo movemos a una columna real.
df_final.reset_index(inplace=True)
df_final.rename(columns={'index': 'sample_id'}, inplace=True)

# 8. GUARDAR EL RESULTADO
# Usamos index=False porque 'sample_id' ya es una columna normal
df_final.to_csv(path_salida, index=False)

print(f"Proceso completado.")
print(f"Dimensiones finales: {df_final.shape[0]} muestras x {df_final.shape[1]} genes/transcritos.")
