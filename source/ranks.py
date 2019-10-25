import pandas as pd 
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

metadata = pd.read_csv("Metadata.csv", index_col = 'Sample' )
strain_occurances = pd.read_csv("Strain_Occurances.csv")
kmer_count = pd.read_csv("kmer_counts.csv", index_col = 'Unnamed: 0')

i = 0 
l = 0

species_occ_dict = {} # key = species, value = # of samples of that species
while l < len(strain_occurances):
	species_occ_dict[strain_occurances.iat[l,0]] = strain_occurances["Number of Occurances"].iat[l]
	l+=1

species_gen_dict = {} # key = genome file name, value = species
while i < len(metadata):
	species_gen_dict[metadata.index.values[i]] = metadata["Species"].iat[i]
	i+=1

cpu_count = cpu_count()
rank_cpu = len(set(species_occ_dict))
sp_match_cpu = cpu_count - rank_cpu - 5

rows = []
kmer_dict = {}
for j in range(0,len(kmer_count)): # for each kmer in the dataframe
	row = kmer_count.iloc[j]
	rows.append(row)
	kmer_dict[row.name] = 0

def sp_match (row): # for a given kmer
	row , species = row
	sp_present = 0 # counter for the amount of times a kmer is present in a sequence of the same species as the function input 
	nsp_present = 0 # counter for the amount of times a kmer is present in a sequence of a different species besides the function input
	kmer = row.name
	for k in range(0,len(row)):	# for every entry in the row
		if row[k] != 0: 
			sample = kmer_count.columns[k]
			sample = sample[0:len(sample)-4]
			if species_gen_dict[sample] == species:
				sp_present +=1
			else: 
				nsp_present+=1 
	return(kmer, sp_present, nsp_present)

def rank(species):
	divisor_sp = species_occ_dict[species] # number of samples with a species matching the input species
	divisor_nsp = len(kmer_count.columns) - divisor_sp # number of samples with a species not matching the input species
	zip_list = zip(rows,[species for i in rows])
	with ProcessPoolExecutor(sp_match_cpu) as ppe:
		for kmer, sp_present, nsp_present in ppe.map(sp_match, zip_list):
			rank = abs((sp_present/divisor_sp)-(nsp_present/divisor_nsp))
			kmer_dict[kmer]=rank
	return(kmer_dict, species)

df = pd.DataFrame(index = kmer_count.index.values , columns = set(species_occ_dict))

with ProcessPoolExecutor(rank_cpu) as ppe:
	for kmer_dict, species in ppe.map(rank, set(species_occ_dict)):
		df[species] = kmer_dict.values()
df.to_csv('Ranks.csv')