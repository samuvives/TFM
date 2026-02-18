# represent the original data vs the predicted reconstructed data
# we have to reconstruct the data
# then we do a graph for each view
# we pass the dataframe to a list and we give it as an axis
# the problem is to use an array. I should check first if the index is the same
# the input is patients (rows) vs features (columns)
# the z matrix is patients (rows) vs factors (columns)
# the w matrix is features (rows) vs factors (columns)
# the reconstruction is patients(rows) vs features (columns)
# the thing is that i run this script and w.index != input.columns, but
# z.index == input.index
# so basically the patients are ordered but the features are not
# para ello z.index == input.index y w.index == input.columns 
# al hacer un shape he visto que las features del input y de los weights no son las mismas
# esto ocurre porque filtra, de tal forma que cuando tu obtienes los weights no estaran
# todas las features originales solo las que no ha reducido el modelo
# listinput = np.log10(np.array(listinput) + 1)
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

INPUTDIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT"
MODELDIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K12"
WEIGHTDIR = os.path.join(MODELDIR, "complete_weights")
OUTPUTDIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/scattersinputvsmodel"
os.makedirs(OUTPUTDIR, exist_ok=True)
viewlist = ["EXPRESSION", "Microbiota", "SV_INS", "VC_11", "Lipidomics", "SV_DEL", "SV_INV", "VC_12", "Metabolomics", "SV_DUP", "SV_TRA"]
continuousviews = ["EXPRESSION", "Microbiota", "Lipidomics", "Metabolomics"]

# input
INPUTFILES = {
    "EXPRESSION": "tpmexpression.tsv",
    "Microbiota": "renamed_microbiota_SOFA_case_tumoral.tsv",
    "SV_INS": "SV_INS_Patient_Gene_Matrix.tsv",
    "VC_11": "MPA_GT_1_1.tsv",
    "Lipidomics": "renamed_lipidomics_data_case.tsv",    
    "SV_DEL": "SV_DEL_Patient_Gene_Matrix.tsv", 
    "SV_INV": "SV_INV_Patient_Gene_Matrix.tsv",  
    "VC_12": "MPA_GT_1_2.tsv",
    "Metabolomics": "renamed_Metabolomics_data_case.tsv",  
    "SV_DUP": "SV_DUP_Patient_Gene_Matrix.tsv",
    "SV_TRA": "SV_TRA_Patient_Gene_Matrix.tsv"
    }

# Z matrix
ZPATH = os.path.join(MODELDIR, "complete_factors_Z_K12.csv")
Z = pd.read_csv(ZPATH, index_col=0)

# weights matrices
WEIGHTFILES = {
    "EXPRESSION": "complete_weights_EXPRESSION_K12.csv",
    "Microbiota": "complete_weights_Microbiota_K12.csv",
    "SV_INS": "complete_weights_SV_INS_K12.csv",
    "VC_11": "complete_weights_VC_11_K12.csv",
    "Lipidomics": "complete_weights_Lipidomics_K12.csv",    
    "SV_DEL": "complete_weights_SV_DEL_K12.csv", 
    "SV_INV": "complete_weights_SV_INV_K12.csv",  
    "VC_12": "complete_weights_VC_12_K12.csv",
    "Metabolomics": "complete_weights_Metabolomics_K12.csv",  
    "SV_DUP": "complete_weights_SV_DUP_K12.csv",
    "SV_TRA": "complete_weights_SV_TRA_K12.csv"
    }

def getinputfile(path):
    inputview = pd.read_csv(os.path.join(INPUTDIR, path), sep="\t", index_col=0)
    return inputview

def getmodeled(path, viewinput):
    W = pd.read_csv(os.path.join(WEIGHTDIR, path), index_col=0)

    print("Dimensions of weight file")
    print(W.shape)
    print("First 5 features in weight file")
    print(W.index[:5])

    print("Dimensions of input file")
    print(viewinput.shape)
    print("First 5 features in input file")
    print(viewinput.columns[:5])

    viewinput = viewinput[W.index]

    print("New dimensions of input file")
    print(viewinput.shape)
    print("First 5 features in new input file")
    print(viewinput.columns[:5])

    if W.index.equals(viewinput.columns):
        W = W.values
        W_t = W.T
        viewmodeled = np.dot(Z.values, W_t)
        return viewmodeled, viewinput
    else:
        return None, None

def scatterinputvsmodeled(view, listinput, listmodeled):
    fig, ax = plt.subplots()
    sns.regplot(x=listinput, y=listmodeled, 
        color = "#99582a",
        scatter_kws={'s': 1, 'alpha': 0.5},
        line_kws={"color": "#6f1d1b"})

    if view in continuousviews:
        ax.set_xlabel("ORIGINAL (log scale)")
    else:
        ax.set_xlabel("ORIGINAL")

    ax.set_ylabel("MODELED")
    ax.set_xscale("log")
    ax.set_title(f"Scatter original vs model, view = {view}")
    savepath = os.path.join(OUTPUTDIR, f"scatteroriginalvsmodel_{view}.png")
    plt.savefig(savepath)
    plt.close()


# MAIN LOOP
for view in viewlist:
    inputpath = INPUTFILES[view]
    viewinput = getinputfile(inputpath)

    print("Dimensions of view input")
    print(viewinput.shape)
    print("First 5 patients in view input:")
    print(viewinput.index[:5])

    print("Dimensions of Z matrix")
    print(Z.shape)
    print("First 5 patients in Z matrix:")
    print(Z.index[:5])
    
    if not viewinput.index.equals(Z.index):
        print("z.index != input.index")
        break

    weightpath = WEIGHTFILES[view]
    viewmodeled, viewinput = getmodeled(weightpath, viewinput)
    if viewmodeled is None:
        print("w.index != input.columns")
        break
    else:
        listinput = viewinput.values.flatten().tolist()
        listmodeled = viewmodeled.flatten().tolist()
        scatterinputvsmodeled(view, listinput, listmodeled)

