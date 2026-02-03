# por cada archivo se filtran las de un genotipo, se aplica el filtro de rare variants y se mete en un dataframe de ese genotipo
# se hace un mpa de ese dataframe
# se repite 3 veces, 1 por cada fenotipo
# a lo mejor hay que poner un identificador del genotipo porque mofa flex no admita dos columnas con el mismo nombre entre variantes
# a lo mejor no tienen un valor por las dos y solo tienen un valor.

import pandas as pd
import os

directorypath = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/newdata/VC"

patient_variants = {"0/1": {}, "1/1": {}, "1/2": {}}

def is_rare_af(af_value):
    """
    Mejorado: Maneja puntos (values missing) y asegura que si
    al menos un valor cumple el criterio, la variante se queda.
    """
    if pd.isna(af_value) or af_value == ".":
        return False

    # Separamos por coma (estándar en VCF para múltiples alelos)
    af_list = str(af_value).split(",")
    for af in af_list:
        try:
            # Si es un punto, float() falla; si es un número, comparamos
            if float(af) < 0.01:
                return True
        except ValueError:
            continue
    return False

for file in os.listdir(directorypath):
    if not file.endswith(".txt"):
        continue

    filepath = os.path.join(directorypath, file)
    patient_name = file.replace("_hg19.txt", "")

    # Cargamos solo las columnas necesarias para ahorrar memoria
    required_cols = ["CHROM", "POS", "REF", "ALT", "AF_nfe", "GEN[*].GT"]
    try:
        df = pd.read_csv(filepath, sep="\t", usecols=required_cols)
    except ValueError as e:
        print(f"Error en {file}: {e}")
        continue

    # 1. Aplicamos el filtro de variantes raras primero
    df_rare = df[df["AF_nfe"].apply(is_rare_af)].copy()

    if df_rare.empty:
        continue

    # 2. Limpieza de genotipo
    df_rare["GEN[*].GT"] = df_rare["GEN[*].GT"].astype(str).str.strip()

    # 3. Crear ID único de forma eficiente
    # Esto evita problemas de tipos de datos al concatenar
    df_rare["VARIANT_ID"] = (
        df_rare["CHROM"].astype(str) + "_" +
        df_rare["POS"].astype(str) + "_" +
        df_rare["REF"].astype(str) + "_" +
        df_rare["ALT"].astype(str)
    )

    # 4. Clasificar por genotipo
    for gt in ["0/1", "1/1", "1/2"]:
        variants_set = set(df_rare.loc[df_rare["GEN[*].GT"] == gt, "VARIANT_ID"])
        if variants_set:
            patient_variants[gt][patient_name] = variants_set

# -----------------------------
# Generación de MPA
# -----------------------------
for gt, data in patient_variants.items():
    if not data:
        print(f"No hay variantes para genotipo {gt}")
        continue

    # Unión de todas las variantes únicas encontradas en este genotipo
    all_variants = sorted(list(set().union(*data.values())))
    all_patients = sorted(data.keys())

    # Crear matriz vacía de ceros
    # Usamos int8 para que el archivo final no pese demasiado
    mpa = pd.DataFrame(0, index=all_patients, columns=all_variants, dtype="int8")

    # Llenar la matriz
    for patient, variants in data.items():
        mpa.loc[patient, list(variants)] = 1

    # Guardar
    output_name = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/VC/MPA_GT_{gt.replace('/', '_')}.csv"
    mpa.to_csv(output_name)
    print(f"Guardado: {output_name} | Pacientes: {len(all_patients)} | Variantes: {len(all_variants)}")
