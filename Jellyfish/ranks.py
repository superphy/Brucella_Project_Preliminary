import pandas as pd 
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

metadata = pd.read_csv("../Metadata.csv", index_col = 'Sample' )
strain_occurances = pd.read_csv("../Strain_Occurances.csv")
kmer_count = pd.read_csv("kmer_counts.csv", index_col = 'Unnamed: 0')
i = 0 
l = 0
a = 0
'''
kmer_count_test = kmer_count.iloc[0:50, 0:30]
kmer_count_test.to_csv("kmer_count_test.csv", index = False)

strains = ['Brucella abortus', 'Brucella melitensis', 'Brucella sp', 'Brucella suis']
sp_dict = {}
for key in strains:
	sp_dict[key] = 0

while a< len(kmer_count.columns):
	#incrementing the correct counter for every occurance of a species name within a group 
	sample = kmer_count.columns[a]
	sample = sample[0:len(sample)-3]
	tax_sp = metadata.loc[sample, 'Species']
	match= re.search('Brucella sp.',tax_sp)
	if match != None: 
		tax_sp_ver = 'Brucella sp'
	if match == None:
		tax_sp_ver = tax_sp
	sp_dict[tax_sp_ver]+=1
	a+=1
print (sp_dict)
'''
species_occ_dict = {} # key = species, value = # of samples of that species
while l < len(strain_occurances):
	species_occ_dict[strain_occurances.iat[l,0]] = strain_occurances["Number of Occurances"].iat[l]
	l+=1

species_gen_dict = {} # key = genome file name, value = species
while i < len(metadata):
	species_gen_dict[metadata.index.values[i]] = metadata["Species"].iat[i]
	i+=1

def rank(species):
	divisor_sp = species_occ_dict[species] # number of samples with a species matching the input species
	divisor_nsp = len(kmer_count.columns) - divisor_sp # number of samples with a species not matching the input species
	out_list = []
	for j in range(0,len(kmer_count)): # for each kmer in the dataframe
		sp_present = 0 # counter for the amount of times a kmer is present in a sequence of the same species as the function input 
		nsp_present = 0 # counter for the amount of times a kmer is present in a sequence of a different species besides the function input
		row = kmer_count.iloc[j] 
		for k in range(0,len(row)):	# for every entry in the row
			if row[k] != 0: 
				sample = kmer_count.columns[k]
				sample = sample[0:len(sample)-3]
				if species_gen_dict[sample] == species:
					sp_present +=1
				else: 
					nsp_present+=1 
		rank = abs((sp_present/divisor_sp)-(nsp_present/divisor_nsp))
		out_list.append(rank)
	return(species, out_list)

df = pd.DataFrame(index = kmer_count.index.values , columns = set(species_occ_dict))

with ProcessPoolExecutor(9) as ppe:
	for species, column in ppe.map(rank, set(species_occ_dict)):
		df[species] = column

df.to_csv('test_rank_div2.csv')