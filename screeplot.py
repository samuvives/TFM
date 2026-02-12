import pandas as pd
import matplotlib.pyplot as plt

# 1. Extraemos el R2 por factor (nos devuelve un diccionario por grupo)
r2_dict = mofaflex_obj.get_r2(total=False, ordered=True)

# 2. Seleccionamos el grupo (si solo tienes uno, suele llamarse 'group1' o el nombre que le dieras)
group_name = list(r2_dict.keys())[0]
df_r2 = r2_dict[group_name]

# 3. Calculamos la varianza total por factor (sumando todas las vistas para ese factor)
# Nota: MOFA reporta fracciones (0.1 = 10%), así que multiplicamos por 100
factor_variance = df_r2.sum(axis=1) * 100
variance_cum = factor_variance.cumsum()

# 4. Creamos la gráfica
fig, ax1 = plt.subplots(figsize=(10, 6))

# Barras para la varianza individual
ax1.bar(factor_variance.index, factor_variance, color='skyblue', label='Varianza por Factor')
ax1.set_ylabel('% Varianza Explicada (Individual)')
ax1.set_xlabel('Factores')

# Línea para la varianza acumulada
ax2 = ax1.twinx()
ax2.plot(factor_variance.index, variance_cum, color='red', marker='o', label='Varianza Acumulada')
ax2.set_ylabel('% Varianza Explicada (Acumulada)')

plt.title(f'Scree Plot - Importancia de Factores ({group_name})')
fig.tight_layout()
plt.show()
