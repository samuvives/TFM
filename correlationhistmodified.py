# study the correlation between factors
# guided and non-guided factors
# heatmap with guided and non-guided
# study the non-guided to help stablish the adequate number of factors
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import stats


def fisher_transform(r):
    """
    Transformación de Fisher para coeficientes de correlación.
    Convierte r (entre -1 y 1) a z-score aproximadamente normal.
    """
    # Manejar valores en los límites (±1) para evitar división por cero
    r = np.clip(r, -0.999999, 0.999999)
    return 0.5 * np.log((1 + r) / (1 - r))


def inverse_fisher_transform(z):
    """
    Transformación inversa de Fisher.
    Convierte z-score de vuelta a coeficiente de correlación.
    """
    return (np.exp(2*z) - 1) / (np.exp(2*z) + 1)


def corrfactorshist(df, namecolumn, hist_title, save_path, apply_fisher=False):
    """
    Crea histograma de valores de correlación.
    
    Args:
        df: DataFrame con los datos
        namecolumn: Nombre de la columna con valores de correlación
        hist_title: Título del histograma
        save_path: Ruta para guardar la imagen
        apply_fisher: Si True, aplica transformación de Fisher
    """
    # Obtener los valores de correlación
    corr_values = df[namecolumn].values
    
    # Aplicar transformación de Fisher si se solicita
    if apply_fisher:
        corr_values = fisher_transform(corr_values)
        # Ajustar título
        hist_title = f"{hist_title} - Fisher Transform"
        xlabel = 'Fisher z-value'
        # Bins específicos para valores transformados (aproximadamente entre -3 y 3)
        bin_width = 0.2
        bins = np.arange(-3, 3 + bin_width, bin_width)
        color = 'lightgreen'
        edgecolor = 'darkgreen'
    else:
        # Datos originales
        xlabel = 'Correlation Coefficient (r)'
        bin_width = 0.05
        bins = np.arange(-1, 1 + bin_width, bin_width)
        color = 'skyblue'
        edgecolor = 'black'
    
    # Crear la figura
    plt.figure(figsize=(12, 8))
    
    # Histograma
    sns.histplot(corr_values, bins=bins, kde=True, 
                 color=color, edgecolor=edgecolor, alpha=0.7)
    
    # Estadísticas descriptivas
    mean_val = np.mean(corr_values)
    std_val = np.std(corr_values)
    skew_val = stats.skew(corr_values)
    kurt_val = stats.kurtosis(corr_values)
    
    # Añadir estadísticas como texto en la gráfica
    stats_text = (f"Mean: {mean_val:.3f}\n"
                  f"Std Dev: {std_val:.3f}\n"
                  f"Skewness: {skew_val:.3f}\n"
                  f"Kurtosis: {kurt_val:.3f}")
    
    # Test de normalidad solo para datos transformados con Fisher
    if apply_fisher and len(corr_values) < 5000:  # Shapiro funciona mejor con n < 5000
        shapiro_stat, shapiro_p = stats.shapiro(corr_values)
        stats_text += f"\nShapiro-Wilk p-value: {shapiro_p:.4f}"
        
        # Interpretación del test
        if shapiro_p > 0.05:
            normality = "NORMAL distribution"
        else:
            normality = "NOT normal distribution"
        
        stats_text += f"\n{normality}"
    
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes,
             fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Añadir línea vertical en la media
    plt.axvline(x=mean_val, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Configuración de la gráfica
    plt.title(hist_title, fontsize=14, fontweight='bold')
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    
    if not apply_fisher:
        # Para datos originales, mostrar marcas específicas
        plt.xticks(np.arange(-1, 1.1, 0.2), rotation=45, ha="right")
    else:
        # Para datos transformados
        plt.xticks(np.arange(-3, 3.1, 0.5), rotation=45, ha="right")
    
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    # Guardar
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {save_path}")
    return corr_values


def create_qq_plot(corr_values, title, save_path):
    """
    Crea un QQ plot para verificar normalidad.
    """
    plt.figure(figsize=(10, 8))
    
    # QQ plot
    stats.probplot(corr_values, dist="norm", plot=plt)
    
    # Añadir línea de referencia y = x
    xlim = plt.xlim()
    ylim = plt.ylim()
    min_val = min(xlim[0], ylim[0])
    max_val = max(xlim[1], ylim[1])
    plt.plot([min_val, max_val], [min_val, max_val], 
             'r--', alpha=0.5, label='y = x (perfect normality)')
    
    plt.title(f'QQ Plot - {title}', fontsize=14, fontweight='bold')
    plt.xlabel('Theoretical Quantiles', fontsize=12)
    plt.ylabel('Sample Quantiles', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def processcorrmatrix(corr_matrix):
    # 1. Limpiar espacios en los nombres para evitar errores visuales
    corr_matrix.columns = corr_matrix.columns.str.replace(' ', '')
    corr_matrix.index = corr_matrix.index.str.replace(' ', '')
    
    # 2. Obtener solo el triángulo superior sin la diagonal (k=1)
    mask = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
    
    # 3. "Apilar" la matriz para pasar de tabla a lista de pares
    df = corr_matrix.where(mask).stack().reset_index()
    df.columns = ['F1', 'F2', 'corrbarvalues']
    
    # 4. Crear el nombre del par para el eje X
    df['corrbarnames'] = df['F1'] + " vs " + df['F2']
    
    return df


klist = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
correlationsumdict = {}
correlationmeandict = {}
fisher_stats = {'sum': {}, 'mean': {}, 'shapiro_p': {}}

print(f"{'='*60}")
print("FISHER TRANSFORMATION ANALYSIS")
print(f"{'='*60}")

for k in klist:
    K_NUMBER = "K" + str(k)
    INPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/MOFAFLEX_FINAL_ANALYSIS/{K_NUMBER}"
    OUTPUT_PATH = f"/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/simpleapproachfinal/postanalysis/histcorrelationanalysis"
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    
    # Crear subdirectorio para gráficas transformadas
    FISHER_OUTPUT_PATH = os.path.join(OUTPUT_PATH, "fisher_transformed")
    os.makedirs(FISHER_OUTPUT_PATH, exist_ok=True)
    
    # Crear subdirectorio para QQ plots
    QQ_OUTPUT_PATH = os.path.join(OUTPUT_PATH, "qq_plots")
    os.makedirs(QQ_OUTPUT_PATH, exist_ok=True)

    print(f"\nProcessing {K_NUMBER}...")

    # Load factors file
    file_z = os.path.join(INPUT_PATH, f"complete_factors_Z_{K_NUMBER}.csv")
    df_Z = pd.read_csv(file_z, index_col=0)

    # Calculamos la correlación de todos los factores
    corr_matrix = df_Z.corr(method='spearman')

    # extraemos aparte solo los no guiados
    corr_matrix_ng = corr_matrix.filter(regex="^Factor")
    corr_matrix_ng = corr_matrix_ng.filter(regex="^Factor", axis=0)

    # procesamos para más visualizaciones
    valuesdf = processcorrmatrix(corr_matrix_ng)
    
    print(f"  Number of correlation pairs: {len(valuesdf)}")
    
    # 1. Histograma original (sin transformación)
    original_title = f"Correlation Distribution - {K_NUMBER}"
    save_path_original = os.path.join(OUTPUT_PATH, f"correlation_hist_original_{K_NUMBER}.png")
    original_values = corrfactorshist(valuesdf, "corrbarvalues", 
                                       original_title, save_path_original,
                                       apply_fisher=False)
    
    # 2. Histograma con transformación de Fisher
    fisher_title = f"Correlation Distribution - {K_NUMBER}"
    save_path_fisher = os.path.join(FISHER_OUTPUT_PATH, f"correlation_hist_fisher_{K_NUMBER}.png")
    
    # Crear histograma con transformación de Fisher
    fisher_values = corrfactorshist(valuesdf, "corrbarvalues",
                                     fisher_title, save_path_fisher,
                                     apply_fisher=True)
    
    # 3. QQ plot para datos transformados con Fisher
    qq_title = f"{K_NUMBER} - Fisher Transformed Values"
    qq_save_path = os.path.join(QQ_OUTPUT_PATH, f"qq_plot_fisher_{K_NUMBER}.png")
    create_qq_plot(fisher_values, qq_title, qq_save_path)
    
    # 4. Guardar estadísticas de Fisher
    fisher_stats['sum'][K_NUMBER] = np.sum(fisher_values ** 2)
    fisher_stats['mean'][K_NUMBER] = np.mean(fisher_values ** 2)
    
    # Test de Shapiro para normalidad
    if len(fisher_values) < 5000:
        _, shapiro_p = stats.shapiro(fisher_values)
        fisher_stats['shapiro_p'][K_NUMBER] = shapiro_p
        normality = "NORMAL" if shapiro_p > 0.05 else "NOT NORMAL"
        print(f"  Shapiro-Wilk test p-value: {shapiro_p:.4f} ({normality})")
    
    # 5. Estadísticas originales (para comparación)
    total_suma = (valuesdf['corrbarvalues'] ** 2).sum()
    correlationsumdict[K_NUMBER] = total_suma
    
    promedio = (valuesdf['corrbarvalues'] ** 2).mean()
    correlationmeandict[K_NUMBER] = promedio
    
    print(f"  Original - Sum of r²: {total_suma:.3f}, Mean of r²: {promedio:.3f}")
    print(f"  Fisher - Sum of z²: {fisher_stats['sum'][K_NUMBER]:.3f}, Mean of z²: {fisher_stats['mean'][K_NUMBER]:.3f}")

# Crear gráficas de resumen
print(f"\n{'='*60}")
print("CREATING SUMMARY PLOTS")
print(f"{'='*60}")

# 1. Gráfica comparativa: Original vs Fisher (suma)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Datos originales
nombres = list(correlationsumdict.keys())
valores_original = list(correlationsumdict.values())
valores_fisher = [fisher_stats['sum'].get(k, 0) for k in nombres]

# Gráfica de suma
ax1.plot(nombres, valores_original, 'o-', color='blue', 
         label='Original (Sum of r²)', linewidth=2, markersize=8)
ax1.plot(nombres, valores_fisher, 's--', color='green', 
         label='Fisher (Sum of z²)', linewidth=2, markersize=8)
ax1.set_title('Comparison: Original vs Fisher Transform (Sum)', fontsize=14, fontweight='bold')
ax1.set_xlabel('Number of factors (K)', fontsize=12)
ax1.set_ylabel('Sum of squared values', fontsize=12)
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.tick_params(axis='x', rotation=45)

# Gráfica de media
valores_original_mean = list(correlationmeandict.values())
valores_fisher_mean = [fisher_stats['mean'].get(k, 0) for k in nombres]

ax2.plot(nombres, valores_original_mean, 'o-', color='blue', 
         label='Original (Mean of r²)', linewidth=2, markersize=8)
ax2.plot(nombres, valores_fisher_mean, 's--', color='green', 
         label='Fisher (Mean of z²)', linewidth=2, markersize=8)
ax2.set_title('Comparison: Original vs Fisher Transform (Mean)', fontsize=14, fontweight='bold')
ax2.set_xlabel('Number of factors (K)', fontsize=12)
ax2.set_ylabel('Mean of squared values', fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.tick_params(axis='x', rotation=45)

plt.tight_layout()
save_path_comparison = os.path.join(OUTPUT_PATH, "original_vs_fisher_comparison.png")
plt.savefig(save_path_comparison, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Comparison plot saved: {save_path_comparison}")

# 2. Gráfica de valores Shapiro-Wilk (normalidad)
if fisher_stats['shapiro_p']:
    plt.figure(figsize=(12, 6))
    
    nombres_shapiro = list(fisher_stats['shapiro_p'].keys())
    valores_shapiro = list(fisher_stats['shapiro_p'].values())
    
    # Crear barras con colores según si son normales o no
    colors = ['green' if p > 0.05 else 'red' for p in valores_shapiro]
    
    bars = plt.bar(nombres_shapiro, valores_shapiro, color=colors, edgecolor='black', alpha=0.7)
    
    # Añadir línea en 0.05 (umbral de significación)
    plt.axhline(y=0.05, color='red', linestyle='--', linewidth=2, alpha=0.5, 
                label='Significance threshold (p=0.05)')
    
    # Añadir etiquetas de valores
    for bar, valor in zip(bars, valores_shapiro):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                 f'{valor:.3f}', ha='center', va='bottom', fontsize=9)
    
    plt.title('Shapiro-Wilk Test p-values for Fisher Transformed Data', 
              fontsize=14, fontweight='bold')
    plt.xlabel('Number of factors (K)', fontsize=12)
    plt.ylabel('p-value', fontsize=12)
    plt.ylim(0, max(valores_shapiro) * 1.2)
    plt.grid(axis='y', alpha=0.3)
    plt.legend()
    plt.xticks(rotation=45)
    
    save_path_shapiro = os.path.join(OUTPUT_PATH, "shapiro_wilk_test_results.png")
    plt.savefig(save_path_shapiro, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Shapiro-Wilk test plot saved: {save_path_shapiro}")

# 3. Crear tabla resumen de estadísticas
print(f"\n{'='*60}")
print("STATISTICAL SUMMARY TABLE")
print(f"{'='*60}")

summary_data = []
for k in klist:
    K_NUMBER = f"K{k}"
    
    # Calcular estadísticas descriptivas básicas
    original_corrs = []
    fisher_zscores = []
    
    # Necesitaríamos volver a cargar los datos para calcular todas las estadísticas
    # Por ahora, usamos las que ya tenemos
    
    row = {
        'K': k,
        'n_pairs': '',  # Se llenará si se quiere calcular
        'Original_Sum_r2': f"{correlationsumdict.get(K_NUMBER, 0):.3f}",
        'Original_Mean_r2': f"{correlationmeandict.get(K_NUMBER, 0):.3f}",
        'Fisher_Sum_z2': f"{fisher_stats['sum'].get(K_NUMBER, 0):.3f}",
        'Fisher_Mean_z2': f"{fisher_stats['mean'].get(K_NUMBER, 0):.3f}",
        'Shapiro_p_value': f"{fisher_stats['shapiro_p'].get(K_NUMBER, 'N/A')}",
        'Is_Normal': 'Yes' if fisher_stats['shapiro_p'].get(K_NUMBER, 0) > 0.05 else 'No'
    }
    
    summary_data.append(row)

summary_df = pd.DataFrame(summary_data)
summary_path = os.path.join(OUTPUT_PATH, "fisher_transformation_summary.csv")
summary_df.to_csv(summary_path, index=False)
print(f"  Summary table saved: {summary_path}")

print(f"\nFirst few rows of summary table:")
print(summary_df.head())

# Mostrar recomendación basada en Fisher
print(f"\n{'='*60}")
print("RECOMMENDATION FOR OPTIMAL NUMBER OF FACTORS")
print(f"{'='*60}")

# Basado en Fisher (queremos menor correlación media)
if fisher_stats['mean']:
    # Buscar K con menor media de z² (menor correlación después de transformación)
    best_k_fisher = min(fisher_stats['mean'], key=fisher_stats['mean'].get)
    best_value_fisher = fisher_stats['mean'][best_k_fisher]
    
    # También considerar normalidad
    best_shapiro = fisher_stats['shapiro_p'].get(best_k_fisher, 0)
    
    print(f"Based on Fisher transformation:")
    print(f"  Optimal K: {best_k_fisher}")
    print(f"  Mean squared z-value: {best_value_fisher:.4f}")
    print(f"  Shapiro-Wilk p-value: {best_shapiro:.4f}")
    
    if best_shapiro > 0.05:
        print(f"  Distribution IS normal after Fisher transform ✓")
    else:
        print(f"  Distribution is NOT normal after Fisher transform")
    
    # Mostrar todos los valores ordenados
    print(f"\nAll K values ranked by Fisher mean squared z-value (lower is better):")
    sorted_k = sorted(fisher_stats['mean'].items(), key=lambda x: x[1])
    for k, value in sorted_k:
        shapiro_val = fisher_stats['shapiro_p'].get(k, 'N/A')
        print(f"  {k}: Mean z² = {value:.4f}, Shapiro p = {shapiro_val}")

print(f"\n{'='*60}")
print("ANALYSIS COMPLETED")
print(f"All plots saved in: {OUTPUT_PATH}")
print(f"{'='*60}")
