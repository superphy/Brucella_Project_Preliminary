import pandas as pd
from Bio import SeqIO

strain_occ = pd.read_csv('../Strain_Occurances.csv', index_col=0)
species = strain_occ.index.values
metadata = pd.read_csv('../Metadata.csv', index_col=0)

def longest(species):
	species_meta = metadata[metadata['Species']==species]
	samples = species_meta.index.values
	longest_counter = 0
	for i in range(0,len(species_meta)):
		sample_len_counter = 0 
		location = '../Approved_Sequences/'+samples[i]+'.fna'
		for record in SeqIO.parse(location, 'fasta'):
			sample_len_counter+=len(record.seq)
		if sample_len_counter > longest_counter:
			longest_counter = sample_len_counter
	return(species, longest_counter)

df = pd.DataFrame(index = species, columns = ['Longest Genome'], data=range(0,len(species)))

for species in species:
	species, longest_counter = longest(species)
	df.at[species, 'Longest Genome' ] = longest_counter

df.to_csv('../Longest_Genomes.csv')