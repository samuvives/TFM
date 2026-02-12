import os
import numpy as np
/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/APPROACHSUBSETS/MOFAFLEX_FINAL_ANALYSIS/K12/complete_weights


def calc_rmse(y_true, y_pred):
    """
    Calcula el Root Mean Square Error entre dos matrices.
    
    Parámetros:
    y_true: Array de NumPy con los datos originales (AnnData.X).
    y_pred: Array de NumPy con la reconstrucción (Z * W.T).
    """
    # 1. Aseguramos que sean arrays de numpy
    y_true = np.asanyarray(y_true)
    y_pred = np.asanyarray(y_pred)
    
    # 2. Calculamos la diferencia al cuadrado
    # Si hay NaNs en tus datos originales, usamos np.nanmean para que no rompa
    squared_errors = (y_true - y_pred) ** 2
    
    # 3. Media de los errores y raíz cuadrada
    mse = np.nanmean(squared_errors)
    rmse = np.sqrt(mse)
    
    return rmse



# real
inputfiles = [f for f in os.listdir(INPUTPATH) if f.endswith(.tsv)]

INPUTPATH =  /gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT
MPA_GT_0_1.tsv                  SV_DUP_Patient_Gene_Matrix.tsv                   renamed_microbiota_SOFA_case_tumoral.tsv
MPA_GT_1_1.tsv                  SV_INS_Patient_Gene_Matrix.tsv                   tpmexpression.tsv
MPA_GT_1_2.tsv                  SV_INV_Patient_Gene_Matrix.tsv    renamed_Metabolomics_data_case.tsv
SV_DEL_Patient_Gene_Matrix.tsv  SV_TRA_Patient_Gene_Matrix.tsv    renamed_lipidomics_data_case.tsv


rmsefiles
views = [EXPRESSION, Microbiota, SV_INS, VC_11, Lipidomics, SV_DEL, SV_INV, VC_12, Metabolomics, SV_DUP, SV_TRA]
# ejecuta para cada view
for view in palabras:

# modeled
WEIGHTSPATH: f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/APPROACHSUBSETS/MOFAFLEX_FINAL_ANALYSIS/K{k}/complete_weights"
complete_weights_EXPRESSION_K12.csv    complete_weights_Microbiota_K12.csv  complete_weights_SV_INS_K12.csv  complete_weights_VC_11_K12.csv
complete_weights_Lipidomics_K12.csv    complete_weights_SV_DEL_K12.csv      complete_weights_SV_INV_K12.csv  complete_weights_VC_12_K12.csv
complete_weights_Metabolomics_K12.csv  complete_weights_SV_DUP_K12.csv      complete_weights_SV_TRA_K12.csv

# busca en cada kfactors file
for k in ks:
    WEIGHTSPATH: f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/APPROACHSUBSETS/MOFAFLEX_FINAL_ANALYSIS/K{k}/complete_weights"
    weightfiles = os.listdir(WEIGHTSPATH)
    for weightfile in weightfiles:
        if palabra in weightfile:
            weightfilepath = os.path.join(WEIGHTSPATH, weightfile)


# creas un rmsefiles por cada view
# cada clave, valor del rmsefiles es un numero de factores
# un diccionario de un diccionario es mejor creo
for key, item in rmsefiles:
    for real, modeled in item:
        # aplicas funcion
        rmseview = calc_rmse(real, modeled)
    rmseviewdict[key] = rmseview

/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUTSUBSETS



results = []
ks = [5, 10, 15, 20, 25] # Los números de factores que quieres probar

for k in ks:
    # 1. Entrenar modelo con k factores
    model = train_sofa(n_factors=k, ...) 
    
    # 2. Predecir para cada vista
    for view_name in ["Metabolomics", "Lipidomics", "Microbiota"]:
        y_true = Xmdata[view_name].X.toarray()
        y_pred = model.predict(site=view_name) # Asumiendo que SOFA permite predecir por sitio
        
        rmse_val = calc_rmse(y_true, y_pred)
        
        results.append({
            'k': k,
            'view': view_name,
            'rmse': rmse_val
        })



import numpy as np
from sofa.utils.utils import calc_rmse

# Supongamos que tienes:
# Z: matriz de factores (n_muestras, n_factores)
# W_metabolomics: matriz de pesos (n_features, n_factores)
# X_original: el array de tu AnnData (n_muestras, n_features)

# 1. Reconstruir la vista específica
# Usamos el producto punto (dot product)
X_pred = np.dot(Z, W_metabolomics.T)

# 2. Calcular el RMSE
# La función comparará celda a celda
error_vista = calc_rmse(X_original, X_pred)

print(f"El error promedio de reconstrucción es: {error_vista}")

# Luego graficas usando Seaborn o Matplotlib
import seaborn as sns
df_res = pd.DataFrame(results)
sns.lineplot(data=df_res, x='k', y='rmse', hue='view', marker='o')
