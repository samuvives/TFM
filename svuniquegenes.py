import pandas as pd
import os

totalgenes= set()

def obtaingenelist(file):
    svfile = pd.read_csv(file, sep="\t") # lee solo las columnas necesarias
    svfile = svfile[svfile["Annotation_mode"] == "split"]
    genes = list(svfile["Gene_name"])
    return genes

path = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/Data/Datos_RNA_Seq_DNA_Guerau/init_15/DNAseq/SV"
for file in [f for f in os.listdir(path) if f.endswith(".tsv")]:
    fullpath = os.path.join(path, file)
    genesinfile = obtaingenelist(fullpath)
    totalgenes.update(genesinfile)

with open("svuniquegenes.txt", "w") as f:
    for gene in totalgenes:
        f.write(f"{gene}\n")


