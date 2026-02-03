def clinicalcorrelationheatmap(Z, gv_names, k, OUTPUT_DIR):
    # create matrix of correlation
    clinical_corr = pd.DataFrame(index=Z.columns, columns=gv_names) # serian los 14 factores como filas y las 2 gv como columnas
    for f_col in Z.columns:
        for g_col in gv_names:
            clinical_corr.loc[f_col, g_col] = Z.loc[:, f_col].corr(gv_df.loc[:, g_col], method='spearman')

    clinical_corr.columns = "GV_" + clinical_corr.columns

    print("clinical_corr head is:")
    print(clinical_corr.head(3))

    # plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(clinical_corr.astype(float), cmap="PiYG", center=0, annot=True, fmt=".2f")
    plt.title(f"Clinical Correlation Heatmap (K:{k})")
    plt.tight_layout()
    savedfilename = os.path.join(OUTPUT_DIR, f"clinical_correlation_heatmap_K{k}.png")
    plt.savefig(savedfilename)
    plt.close()

    print("Clinical Correlation Heatmap done")
    print(f"Result saved in path: {savedfilename}")

if __name__ == "main":
    # graph of clinical correlation (factors vs guiding variables)
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    k = "12"
    INPUT_BASE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/"
    INPUT_DIR = os.path.join(INPUT_BASE, f"K{k})
    GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"
    OUTPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis"

    # load factors
    factors_file = os.path.join(INPUT_DIR, f"complete_factors_Z_K{k}.csv")
    Z = pd.read_csv(factors_file, index_col=0)

    print("Z head is:")
    print(Z.head(3))

    # load gv
    gv_df = pd.read_csv(GV_FILE, index_col=0)
    gv_names = gv_df.columns.tolist()

    print("gv file head is:")
    print(gv_df.head(3))
    clinicalcorrelationheatmap(Z, gv_names)

