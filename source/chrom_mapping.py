import pandas as pd 
import numpy as np
import ahocorasick
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

metadata = pd.read_csv('../Metadata.csv')
tkmers = pd.read_csv('../Ranks.csv', index_col = 0)
strain_occurances = pd.read_csv("../Strain_Occurances.csv", index_col = 0)
#species = strain_occurances.index.values 
species = ['Brucella melitensis']
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a' }
meta_ref = metadata[metadata['Number of Contigs'] == 2] # only complete sequences (a.k.a 2 contigs)
meta_ref = meta_ref.drop_duplicates(subset = 'Species', keep='first') # one of each species

c1_index, c1_rank, c2_index, c2_rank = {}, {}, {}, {}
for kmer in tkmers.index: # for each kmer in the dataframe
		c1_index[kmer]=0
		c1_rank[kmer]=0
		c2_index[kmer]=0
		c2_rank[kmer]=0

#Choosing and locating a refrence genome 
def refrence_location(species):
	ref_row = meta_ref[meta_ref['Species']==species]
	ref_file = ref_row['Sample'].to_list()
	ref_file_location = '../Approved_Sequences/'+ref_file[0]+'.fna'
	return(ref_file_location)

#Creates 2 haystacks from the refrence genome, one for each chromosome
def haystacks(species):
	location = refrence_location(species)
	sequence = []
	with open(location) as handle:
		for title, seq in SimpleFastaParser(handle):
			sequence.append(seq)
	chromosome1 = sequence[0]
	chromosome2 = sequence[1]
	return(chromosome1, chromosome2)

#Creates the needles, which are the top X kmers identified for each species
def needles(species):
	kmers = tkmers.index.values
	rank = tkmers[species+' RANK']
	aho = ahocorasick.Automaton()
	for i in range(0,len(kmers)):
		aho.add_word(kmers[i],(kmers[i], rank[kmers[i]]))
	aho.make_automaton()
	return(aho)

#Runs the ahocorasick algorithm to find the location of the kmers (needles) in the refrence genome (haystack) 
def find_needles(species):
	aho = needles(species)
	chromosome1, chromosome2 = haystacks(species)
	c1_index, c1_rank, c2_index, c2_rank = [], [], [], []
	for end_index, (kmer, rank) in aho.iter(chromosome1):
		c1_index.append(end_index)
		c1_rank.append(rank)
	for end_index, (kmer, rank) in aho.iter(chromosome2):
		c2_index.append(end_index)
		c2_rank.append(rank)
	return(species, c1_index, c1_rank, c2_index, c2_rank)


with ProcessPoolExecutor(len(species)) as ppe:
	for species, c1_index, c1_rank, c2_index, c2_rank in ppe.map(find_needles, species):
		fig, (ax1, ax2) = plt.subplots(1,2, sharey = True)
		fig.suptitle(species)

		c1_dict, c2_dict = {}, {}
		c1_length = len(max(c1_index)-min(c1_index))
		c1_index_padded = range(min(c1_index), max(c1_index))
		c1_rank_padded = np.zeros(c1_length)

		for i in range(0,len(c1_rank)):
			index = c1_index[i]
			c1_rank_padded[index]=c1_rank[i]

		c2_length = len(max(c2_index)-min(c2_index))
		c2_index_padded = range(min(c2_index), max(c2_index))
		c2_rank_padded = np.zeros(c2_length)

		for i in range(0,len(c2_rank)):
			index = c2_index[i]
			c2_rank_padded[index]=c2_rank[i]
'''
		for i in range(min(c1_index), max(c1_index), 1):
			c1_dict[i]=0
		for i in range(min(c2_index), max(c2_index), 1):
			c2_dict[i]=0
		for i in range(0,len(c1_rank)):
			c1_dict[c1_index[i]]=c1_rank[i]
		for i in range(0,len(c2_rank)):
			c2_dict[c2_index[i]]=c2_rank[i]
	
		df_c1_rank = pd.DataFrame(c1_dict.values())
		roll_c1 = df_c1_rank.rolling(window=1000).mean()
		df_c2_rank = pd.DataFrame(c2_dict.values())
		roll_c2 = df_c2_rank.rolling(window=1000).mean()
		
		c1_keys = []
		for key in c1_dict.keys():
			c1_keys.append(key)
		
		c2_keys = []
		for key in c2_dict.keys():
			c2_keys.append(key)
'''		
		df_c1_rank = pd.DataFrame(c1_rank_padded)
		roll_c1 = df_c1_rank.rolling(window=1000).mean()
		df_c2_rank = pd.DataFrame(c2_rank_padded)
		roll_c2 = df_c2_rank.rolling(window=1000).mean()
		
		ax1.plot(c1_index_padded, c1_rank_padded, color= '#cbc9cc', linewidth=0.25, ls = ':')
		ax1.plot(c1_index_padded, roll_c1, color=col_dict[species], linewidth=1.5)
		#bins_c1 = range(0,max(c1_index),1000)
		#ax1.hist(c1_index, bins = bins_c1,color=col_dict[species])
		#ax1.scatter(c1_index, c1_rank, color=col_dict[species], marker='.')
		#ax1.set_ylim([-1,1])
		ax1.set_title('Chromosome 1')
		
		
		ax2.plot(c2_index_padded, c2_rank_padded, color='#cbc9cc', linewidth=0.25, ls = ':')
		ax2.plot(c2_index_padded, roll_c2, color=col_dict[species], linewidth=1.5)
		#bins_c2 = range(0,max(c2_index),1000)
		#ax2.hist(c2_index, bins = bins_c2,color=col_dict[species])
		#ax2.scatter(c2_index, c2_rank, color=col_dict[species], marker='.')
		#ax2.set_ylim([-1,1])
		ax2.set_title('Chromosome 2')
		
		plt.show()