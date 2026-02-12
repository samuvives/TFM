
    # ===================================================
    # 3.14 CLINIC CORRELATION (FACTORS VS GV)
    # ===================================================
    print("Calculating clinic correlation")

    from clinicalcorrelationmodified import clinicalcorrelationheatmap
    INPUT_BASE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/"
    INPUT_DIR = os.path.join(INPUT_BASE, k_number)
    GV_FILE = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/MOFAINPUT/matrizGV4_mapeado.csv"
    OUTPUT_DIR = "/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis"

    # load factors
    factors_file = os.path.join(INPUT_DIR, f"complete_factors_Z_{k_number}.csv")
    Z = pd.read_csv(factors_file, index_col=0)

    print("Z head is:")
    print(Z.head(3))

    # load gv
    gv_df = pd.read_csv(GV_FILE, index_col=0)
    gv_names = gv_df.columns.tolist()

    print("gv file head is:")
    print(gv_df.head(3))
    clinicalcorrelationheatmap(Z, gv_names, k, OUTPUT_DIR)
