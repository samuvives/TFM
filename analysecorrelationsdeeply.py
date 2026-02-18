# we compute the number of times the correlation between factors is significant
# pvalue < 0.05
# correlation = rho
# how do you apply the spearmanr to every combination of factors
# count the number of factors 
# lo que haces es que cuentas el número de factores que tienes (columnas en tu dataframe)
# creas dataframe para contener los valores y otro para contener los p values

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

approach = "simpleapproachfinal"

PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/{approach}/"
INPUTDIRPATH = os.path.join(PATH, "MOFAFLEX_FINAL_ANALYSIS")
OUTPUTPATH = os.path.join(PATH, "postanalysis/numsignificantcorrelationsperrun.png")
factordirs = [f for f in os.listdir(INPUTDIRPATH) if f.startswith("K")]


def getmatrices(Z):
    n = len(Z.columns)
    correlationmatrix = pd.DataFrame(np.zeros((n, n)), index=Z.columns, columns=Z.columns)
    pvaluematrix = pd.DataFrame(np.zeros((n, n)), index=Z.columns, columns=Z.columns)

    for i in range(n):
        for j in range(n):
            correlation, pvalue = spearmanr(Z.iloc[:, i], Z.iloc[:, j])
            correlationmatrix.iloc[i, j] = correlation
            pvaluematrix.iloc[i, j] = pvalue
    return correlationmatrix, pvaluematrix

def matrixtotable(df, colname):
    mask = np.triu(np.ones_like(df, dtype=bool), k=1)
    df = df.where(mask)
    df = df.stack().reset_index()
    df.columns = ["FactorA", "FactorB", colname]
    df["FactorvsFactor"] = df["FactorA"] + "vs" + df["FactorB"]
    colstodrop = ["FactorA", "FactorB"]
    df = df.drop(columns=colstodrop)
    return df

runsdict = {}

for factordir in factordirs:
    INPUTPATH = os.path.join(PATH, f"MOFAFLEX_FINAL_ANALYSIS/{factordir}/complete_factors_Z_{factordir}.csv")
    Z = pd.read_csv(INPUTPATH, index_col=0)

    correlationmatrix, pvaluematrix = getmatrices(Z)

    correlationtable = matrixtotable(correlationmatrix, "Correlation")
    pvaluetable = matrixtotable(pvaluematrix, "p-value")
    factorstable = pd.merge(correlationtable, pvaluetable, on="FactorvsFactor")
    print(factorstable)

    numbersignificantcorrelations = (factorstable["p-value"] < 0.05).sum()
    runsdict[factordir] = numbersignificantcorrelations

dfruns = pd.DataFrame.from_dict(runsdict, orient="index").reset_index()
dfruns.columns = ["Runwithxfactors", "Numberofsignificantcorrelations"]
dfruns["Runnumber"] = dfruns["Runwithxfactors"].str.replace("K", "").astype(int)
dfruns = dfruns.sort_values(by="Runnumber")

fig, ax = plt.subplots()
ax.scatter(dfruns["Runwithxfactors"], dfruns["Numberofsignificantcorrelations"], color= "#fb8500")
ax.plot(dfruns["Runwithxfactors"], dfruns["Numberofsignificantcorrelations"], color= "#fb8500")
ax.set_xlabel("Run with x factors")
ax.set_ylabel("Number of significant correlations")
ax.set_title("Number of significant correlations across the runs")
plt.savefig(OUTPUTPATH)
plt.close()
