import pandas as pd 
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

ranks = pd.read_csv('Ranks.csv', index_col = 0)
species = ranks.columns.values 
sample_size = 100

kmer_dict = {}
for kmer in ranks.index: # for each kmer in the dataframe
	kmer_dict[kmer] = np.inf

def sort(species):
	species_df = abs(ranks[species])
	species_df = species_df.sort_values(ascending=False)
	species_df = species_df.iloc[0:sample_size]
	for i in range(0,len(species_df)):
		kmer = species_df.index[i]
		value = species_df[i]
		kmer_dict[kmer]=value
	return (kmer_dict, species)

df = pd.DataFrame(index = kmer_dict, columns = species)

with ProcessPoolExecutor(cpu_count()-10) as ppe:
	for kmer_dict, species in ppe.map(sort, species):
		df[species] = kmer_dict.values()

df = df.dropna(axis='index', how='all')
df.to_csv('Top_kmers.csv')


