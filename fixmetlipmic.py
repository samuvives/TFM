# convert from csv to tsv and substitute the spaces in the colnames with underscores
import pandas as pd
import os

directory = '/gpfs/projects/bsc20/bsc236340/Project_IDIBAPS/NEWMETLIPMIC'

files = os.listdir(directory)

for name_file in files:
    if name_file.endswith('.csv'):
        input_path = os.path.join(directory, name_file)
        try:
            df = pd.read_csv(input_path)
            
            df.columns = [col.strip().replace(' ', '_').replace('"', '') for col in df.columns]
            
            name_tsv = name_file.rsplit('.', 1)[0] + '.tsv'
            output_path = os.path.join(directory, name_tsv)
            df.to_csv(output_path, sep='\t', index=False)
            print(f"Saved in: {output_path}")
            
        except Exception as e:
            print(f"Error en {name_file}: {e}")

print("\n Process finished")
