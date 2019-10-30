import pandas as pd 
import numpy as np
import ahocorasick
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

metadata = pd.read_csv('../Metadata.csv')
tkmers = pd.read_csv('../Top_kmers.csv', index_col = 0)

species = tkmers.columns.values[0:len(tkmers.columns.values)]
meta_ref = metadata[metadata['Number of Contigs'] == 2] # only complete sequences (a.k.a 2 contigs)
meta_ref = meta_ref.drop_duplicates(subset = 'Species', keep='first') # one of each species
chroms = [' c1', ' c2']
columns = [y+x for y in species for x in chroms]

chrom1_dict = {}
chrom2_dict = {}
for kmer in tkmers.index:
	chrom1_dict[kmer] = np.inf
	chrom2_dict[kmer] = np.inf

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
	kmers = tkmers[species].dropna(axis='index', how='all').index.values
	aho = ahocorasick.Automaton()
	for i in range(0,len(kmers)):
		aho.add_word(kmers[i],kmers[i])
	aho.make_automaton()
	return(aho)

#Runs the ahocorasick algorithm to find the location of the kmers (needles) in the refrence genome (haystack) 
def find_needles(species):
	aho = needles(species)
	chromosome1, chromosome2 = haystacks(species)
	for end_index, kmer in aho.iter(chromosome1):
		chrom1_dict[kmer] = end_index
	for end_index, kmer in aho.iter(chromosome2):
		chrom2_dict[kmer] = end_index
	return(species, chrom1_dict, chrom2_dict)

df = pd.DataFrame(columns = columns, index = tkmers.index.values)

with ProcessPoolExecutor(len(species)) as ppe:
	for species, chrom1_dict, chrom2_dict in ppe.map(find_needles, species):
		df[species+' c1'] = chrom1_dict.values()
		df[species+' c2'] = chrom2_dict.values()

df.to_csv('Kmers_Locations.csv')
print(df)