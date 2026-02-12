# check if the factors are normally distributed
import os
import numpy as np
import pandas as pd
from scipy import stats
filepath = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K12/complete_factors_Z_K12.csv"
Zmatrix = pd.read_csv(filepath, index_col=0)
def checknormal(factor, factorname):
    res = stats.normaltest(factor)
    pvaluefactor = res.pvalue
    if pvaluefactor > 0.05:
        print(f"{factorname} follows a normal distribution")
    else:
        print(f"{factorname} does not follow a normal distribution")

for factorname in Zmatrix:
    factor = Zmatrix[factorname]
    checknormal(factor, factorname)
