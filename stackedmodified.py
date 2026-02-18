# cada archivo una view, cada fila una feature, cada columna un factor
import os
import pandas as pd
import matplotlib.pyplot as plt


approach = "OBTAININGELBO"
PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/{approach}"
INPUTDIRPATH = os.path.join(PATH, "MOFAFLEX_FINAL_ANALYSIS")
LISTDIRS = [f for f in os.listdir(INPUTDIRPATH) if os.path.isdir(os.path.join(INPUTDIRPATH, f))]
OUTPUTDIR = os.path.join(PATH, "postanalysis")
os.makedirs(OUTPUTDIR, exist_ok=True)

def stackedfactorsactivity(df_r2, KDIR, OUTPUTPATH):
    df_r2 = df_r2.copy()
    # barplot stacked
    df_r2['Total_Var'] = df_r2.sum(axis=1)

    # 2. Ordenar de mayor a menor y eliminar la columna temporal para el gráfico
    df_r2 = df_r2.sort_values(by='Total_Var', ascending=False).drop(columns=['Total_Var'])
    df_r2_nonguided = df_r2[df_r2.index.str.startswith("Factor")]


    ax = df_r2.plot(kind='bar',
                  stacked=True,
                  figsize=(12, 7),
                  colormap="tab20",
                  edgecolor='white',
                  linewidth=0.5)

    plt.title(f"Factor composition per view. {KDIR}", fontsize=16, pad=20)
    plt.ylabel("Activity", fontsize=12)
    plt.xlabel("Factors", fontsize=12)
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")

    plt.legend(title="Omic views",
               bbox_to_anchor=(1.05, 1),
               loc='upper left',
               fontsize=10,
               frameon=False)

    plt.tight_layout()
    plt.savefig(OUTPUTPATH, dpi=300, bbox_inches='tight')
    plt.close()

def stackedfactorsactivitywiththr(df_r2, KDIR, OUTPUTPATHWITHTHR):
    df_r2 = df_r2.copy()

    # barplot stacked
    df_r2['Total_Var'] = df_r2.sum(axis=1)

    # 2. Ordenar de mayor a menor y eliminar la columna temporal para el gráfico
    df_r2 = df_r2.sort_values(by='Total_Var', ascending=False).drop(columns=['Total_Var'])
    df_r2_nonguided = df_r2[df_r2.index.str.startswith("Factor")]

    max_height = df_r2_nonguided.sum(axis=1).max()
    threshold = max_height * 0.01
    threshold2 = max_height * 0.02

    ax = df_r2.plot(kind='bar',
                  stacked=True,
                  figsize=(12, 7),
                  colormap="tab20",
                  edgecolor='white',
                  linewidth=0.5)

    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1, label=f'Threshold (1%)')
    ax.axhline(y=threshold2, color='green', linestyle='--', linewidth=1, label=f'Threshold (2%)')
    plt.title(f"Factor composition per view. {KDIR}", fontsize=16, pad=20)
    plt.ylabel("Activity", fontsize=12)
    plt.xlabel("Factors", fontsize=12)
    plt.xticks(rotation=45, ha="right", rotation_mode="anchor")

    plt.legend(title="Omic views",
               bbox_to_anchor=(1.05, 1),
               loc='upper left',
               fontsize=10,
               frameon=False)

    plt.tight_layout()
    plt.savefig(OUTPUTPATHWITHTHR, dpi=300, bbox_inches='tight')
    plt.close()

for KDIR in LISTDIRS:
    fullpath = os.path.join(INPUTDIRPATH, KDIR)
    fullpath = os.path.join(fullpath, "varianza_explicada_por_factor.csv")
    df_r2 = pd.read_csv(fullpath, index_col=0)

    OUTPUTPATH = os.path.join(OUTPUTDIR, f"factor_composition_stacked_{KDIR}.png")
    stackedfactorsactivity(df_r2, KDIR, OUTPUTPATH)

    OUTPUTPATHWITHTHR = os.path.join(OUTPUTDIR, f"factor_composition_stacked_thr_{KDIR}.png")
    stackedfactorsactivitywiththr(df_r2, KDIR, OUTPUTPATHWITHTHR)




