import pandas as pd
import Bio 

strain_occ = pd.read_csv('../Strain_Occurances.csv', index_col=0)
species = strain_occ.index.values
metadata = pd.read_csv('../Metadata.csv', index_col=0)

def longest(species):
	species_meta = metadata[metadata['Species']==species]
	samples = species_meta.index.values
	for i in range(0,len(species_meta)):
		print(samples[i])

longest('Brucella suis')