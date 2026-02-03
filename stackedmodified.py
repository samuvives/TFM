# como se obtuvieron los weights
# cada archivo una view, cada fila una feature, cada columna un factor

def processweightfiles(WEIGHTS_DIR):
    view_activity_dict = {}

    # search files
    # Rows = Genes/Variables, Columns = Factors
    for file in os.listdir(WEIGHTS_DIR):

        # get view name
        viewname = file.replace(".csv", "")
        viewname = viewname.split("_")[2]

        # read file
        W = pd.read_csv(file)
        Wseries = W**2.sum(axis=0)
        view_activity_dict[viewname: Wseries]

    df_activity_all = pd.DataFrame(view_activity_dict)

    # sumas todas las columnas y te quedas con un valor global por factor
    total_activity = (df_activity_all
        .sum(axis=1)
        .sort_values(ascending=False))

    df_plot = df_activity_all.loc[total_activity.index]

    return df_plot

def stackedgraph(df_plot, k):

    # ==========================================
    # 3. GENERACIÓN DEL GRÁFICO STACKED
    # ==========================================
    ax = df_plot.plot(kind='bar', 
                      stacked=True, 
                      figsize=(12, 7), 
                      colormap="tab20", 
                      edgecolor='white', 
                      linewidth=0.5)

    plt.title(f"Factor composition per view. {K_FOLDER}", fontsize=16, pad=20)
    plt.ylabel("Sum of squared weights (Activity)", fontsize=12)
    plt.xlabel("Factores", fontsize=12)
    plt.xticks(rotation=45)

    plt.legend(title="Omic views)", 
               bbox_to_anchor=(1.05, 1), 
               loc='upper left', 
               fontsize=10, 
               frameon=False)

    plt.tight_layout()

    save_path = os.path.join(OUTPUT_PATH, f"factor_composition_stacked_{K_FOLDER}.png")
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    print(f"--- GRÁFICO GENERADO ---")
    print(f"Archivo guardado en: {save_path}")

if __name__ == "__main__":

    import pandas as pd
    import matplotlib.pyplot as plt
    import os

    # parameters
    k = 30
    K_FOLDER = "K" + str(k)

    WEIGHTS_DIR = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_FOLDER}/complete_weights/"
    OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/{K_FOLDER}"
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    df_plot = processweightfiles(WEIGHTS_DIR)
    stackedgraph(df_plot, k)

