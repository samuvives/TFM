import os
import pandas as pd
path = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/K12/complete_weights"
weightslist = [f for f in os.listdir(path) if f.startswith("complete")]
for weightfile in weightslist:
    weighfile = os.path.join(path, weighfile)
    weightdata = pd.read_csv(weightfile)
    print(weightfile)
    weightdata.describe()
    print("----")
