"""a ver el script ha de:
crear un data frame vacio por cada tipo de SV
leer los dataframes de cada tipo
leer el excel con las correspondencias de gen a genename /home/vant/Escritorio/TFM/datostfm/trabajoSamuel/SV/SVgenetypes.xls
leer el archivo, solo las columnas necesarias, es decir, tipo de fila, tipo de sv y genename 
quedarse solo con las que son split, que sean unicas
añadir una columna con la lista de genes unicos a términos hugo de biomart, y dejar los que no tienen correspondencia.
usar la columna de identificador de paciente y de genes para crear la matriz
guardar la matriz en el dataframe correspondiente y generar mediante ese dataframe un archivo"""
# la columna de tipo de fila es Annotation_name
# la columna de tipo de SV es SV_type

import pandas as pd
import os
import openpyxl

# 1. Definir la ruta y obtener lista de archivos
ruta_carpeta = 'tu/ruta/de/archivos/'
archivos_csv = glob.glob(os.path.join(ruta_carpeta, "*.csv"))

lista_df = []

for archivo in archivos_csv:
    # Extraemos solo el nombre del archivo (sin la ruta completa) para la etiqueta
    nombre_base = os.path.basename(archivo)

    # 2. Leemos solo la columna necesaria (usa el nombre real de tu columna)
    df_temp = pd.read_csv(archivo, usecols=['tu_columna_de_datos', ""])

    # 3. Creamos la columna de trazabilidad con el nombre del archivo
    df_temp['archivo_origen'] = nombre_base

    # Añadimos a la lista para luego unir todo
    lista_df.append(df_temp)

# 4. Consolidamos todo en un solo DataFrame final
df_final = pd.concat(lista_df, ignore_index=True)

print(df_final.head())


# la movida es que tengo que extraer del archivo, es una lista por cada paciente de los genes que tienen
# esa lista va a un dataframe de pandas, el paciente es la key del diccionario con el que se crea
def insertgenes(svfile):
    """Pone la correspondencia con la columna Gene_name
    y si no hay correspondencia deja la que está.
    Devuelve el dataframe resultante"""
    # Unimos los dataframes (Left Join)
    df_resultado = pd.merge(df1, df2, left_on='ID', right_on='Original', how='left')

    # Si no hubo coincidencia (NaN), usamos el valor de la columna original
    df_resultado['Nueva_Columna'] = df_resultado['Nuevo'].fillna(df_resultado['ID'])

    # (Opcional) Eliminamos las columnas extra que trajo el merge
    df_resultado = df_resultado.drop(columns=['Original', 'Nuevo'])


def creatempa(df_final):
    df_final['presencia'] = 1
    df_matriz = df_final.pivot(
        index='archivo_origen',
        columns='tu_columna_de_datos',
        values='presencia'
    )
    df_matriz = df_matriz.fillna(0).astype(int)

print(df_matriz)

path = "gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/finalsubsets"
genelist = pd.read_excel(os.path.join(path, "svgenes.xlsx"))
genesperpatient = {}

for file in [f for f in os.listdir(path) if f.endswith(".tsv")]:
    base_name = os.path.basename(archivo)
    fullpath = os.path.join(path, file)
    svfile = pd.read_csv(fullpath, sep="\t", usecols=["", ""])
    svfile = insertgenes(svfile)
    df_temp = sv_file[["officialname"]]
    df_temp['archivo_origen'] = base_name # pones nombre de archivo en todas las filas del dataframe
    lista_df.append(df_temp) # lista de dataframes
df_final = pd.concat(lista_df, ignore_index=True)

