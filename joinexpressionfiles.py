import pandas as pd

# Cargar los archivos indicando que el separador es un tabulador (\t)
df1 = pd.read_csv('archivo1.tsv', sep='\t')
df2 = pd.read_csv('archivo2.tsv', sep='\t')

# Unir lateralmente (axis=1)
resultado = pd.concat([df1, df2], axis=1)
resultado = resultado.T

# Guardar el resultado
resultado.to_csv('unido.tsv', sep='\t', index=False)
