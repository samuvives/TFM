import os
import pandas as pd
PATH = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/MPA_GT_0_1.tsv"
data = pd.read_csv(PATH, sep="\t", index_col=0)
