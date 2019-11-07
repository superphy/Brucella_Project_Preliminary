import pandas as pd 
import re, sys
import scipy.stats as stats
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

if __name__ == '__main__':
	
	metadata = pd.read_csv("../Metadata.csv", index_col = 'Sample' )
	strain_occurances = pd.read_csv("../Strain_Occurances.csv")
	kmer_count = pd.read_csv("../kmer_count_test.csv", index_col = 'Unnamed: 0')

	species_occ_dict = {} # key = species, value = # of samples of that species
	for i in range(0,len(strain_occurances)):
		species_occ_dict[strain_occurances.iat[i,0]] = strain_occurances["Number of Occurances"].iat[i]

	species_gen_dict = {} # key = genome file name, value = species
	for j in range(0,len(metadata)):
		species_gen_dict[metadata.index.values[j]] = metadata["Species"].iat[j]

	a = 0 	
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

	rows = []
	rank_dict = {}
	p_dict = {}
	for k in range(0,len(kmer_count)): # for each kmer in the dataframe
		row = kmer_count.iloc[k]
		rows.append(row)
		rank_dict[row.name] = 0
		p_dict[row.name] = 0

	cpu_count = cpu_count()
	rank_cpu = len(set(sp_dict))
	sp_match_cpu = cpu_count - rank_cpu - 5

	def sp_match(row): # for a given kmer
		row , species = row
		sp_present = 0 # counter for the amount of times a kmer is present in a sequence of the same species as the function input 
		nsp_present = 0 # counter for the amount of times a kmer is present in a sequence of a different species besides the function input
		kmer = row.name
		for l in range(0,len(row)):	# for every entry in the row
			if row[l] != 0: 
				sample = kmer_count.columns[l]
				sample = sample[0:len(sample)-3]
				if species_gen_dict[sample] == species:
					sp_present +=1
				else: 
					nsp_present+=1 
		return(kmer, sp_present, nsp_present)

	def rank(species):
		for key in rank_dict: 
			rank_dict[key] = 0
		for key in p_dict: 
			p_dict[key] = 0
		divisor_sp = sp_dict[species] # number of samples with a species matching the input species
		divisor_nsp = len(kmer_count.columns) - divisor_sp # number of samples with a species not matching the input species
		zip_list = zip(rows,[species for i in rows])
		with ProcessPoolExecutor(sp_match_cpu) as ppe:
			for kmer, sp_present, nsp_present in ppe.map(sp_match, zip_list):
				rank = (sp_present/divisor_sp)-(nsp_present/divisor_nsp)
				odds_ratio, p_value = stats.fisher_exact([[sp_present, nsp_present],[divisor_sp, divisor_nsp]])
				p_dict[kmer]=p_value
				rank_dict[kmer]=rank
		return(rank_dict,p_dict, species)

	df = pd.DataFrame(index = kmer_count.index.values)# , columns = set(sp_dict))

	with ProcessPoolExecutor(rank_cpu) as ppe:
		for rank_dict, p_dict, species in ppe.map(rank, set(sp_dict)):
			df[species+' RANK'] = rank_dict.values()
			df[species+' P_VAL'] = p_dict.values()
	
	print(df)
	df.to_csv('Ranks_test.csv')