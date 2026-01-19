# por cada archivo se filtran las de un genotipo, se aplica el filtro de rare variants y se mete en un dataframe de ese genotipo
# se hace un mpa de ese dataframe
# se repite 3 veces, 1 por cada genotype
# a lo mejor hay que poner un identificador del genotipo porque mofa flex no admita dos columnas con el mismo nombre entre variantes
# a lo mejor no tienen un valor por las dos y solo tienen un valor.

import pandas as pd
import os

directorypath = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/newdata/VC"

patient_variants = {"0/1": {}, "1/1": {}, "1/2": {}}

def is_rare_af(af_value):
    """
    Filters rare variants. Manage dots (values missing) and 
    if any value comply to the criteria the variant stays
    """
    if pd.isna(af_value) or af_value == ".":
        return False

    # divide by comma when multiple altered alleles
    af_list = str(af_value).split(",")
    for af in af_list:
        try:
            # float() fails if its a dot; compares if its a number
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

    # Load only the required columns
    required_cols = ["CHROM", "POS", "REF", "ALT", "AF_nfe", "GEN[*].GT"]
    try:
        df = pd.read_csv(filepath, sep="\t", usecols=required_cols)
    except ValueError as e:
        print(f"Error in {file}: {e}")
        continue

    # Apply rare variants filter
    df_rare = df[df["AF_nfe"].apply(is_rare_af)].copy()

    if df_rare.empty:
        continue

    # 2. Limpieza de genotipo
    df_rare["GEN[*].GT"] = df_rare["GEN[*].GT"].astype(str).str.strip()

    # Create ID
    df_rare["VARIANT_ID"] = (
        df_rare["CHROM"].astype(str) + "_" +
        df_rare["POS"].astype(str) + "_" +
        df_rare["REF"].astype(str) + "_" +
        df_rare["ALT"].astype(str)
    )

    # Classify per genotype
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

    # union of all variants in a genotype
    all_variants = sorted(list(set().union(*data.values())))
    all_patients = sorted(data.keys())

    # Empty matrix of zeros
    mpa = pd.DataFrame(0, index=all_patients, columns=all_variants, dtype="int8")

    # Fill the matrix
    for patient, variants in data.items():
        mpa.loc[patient, list(variants)] = 1

    # Save
    output_name = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/VC/MPA_GT_{gt.replace('/', '_')}.csv"
    mpa.to_csv(output_name)
    print(f"Guardado: {output_name} | Pacientes: {len(all_patients)} | Variantes: {len(all_variants)}")
