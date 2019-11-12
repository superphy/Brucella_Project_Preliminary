import pandas as pd
import numpy as np

og_meta = pd.read_csv('../Metadata.csv', index_col = 0)
strain_occ = pd.read_csv('../Strain_Occurances.csv', index_col = 0)
kcounts = pd.read_csv('../kmer_counts.csv')

index = []
for name in og_meta.index:
	index.append(name+'.fna')

kover_meta = pd.DataFrame(index = index, columns = strain_occ.index, data=np.zeros((len(index),len(strain_occ.index))) )

#Generates a large metadata folder with each column representing a species
for i in range(0, len(og_meta.index)):
	file = index[i]
	species = og_meta['Species'].iloc[i]
	kover_meta.loc[file,species]=1

#Parses the large metadata folder to be one file per species
for i in range(0, len(strain_occ.index)):
	strain = strain_occ.index[i]
	kover_single_meta = kover_meta[strain]
	under_strain = strain.replace(' ', '_')
	kover_single_meta.to_csv('../Kover_Data/Kover_'+under_strain+'_Metadata.tsv', sep='\t', header=False)

#Created the kmer matrix file 
kcounts.rename(columns={'Unnamed: 0':'kmers'}, inplace=True )
kcounts.astype(bool)
kcounts.to_csv('../Kover_Data/Kmer_Matrix.tsv', sep='\t', index=False)
