import pandas as pd 
import numpy as np
import ahocorasick
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

metadata = pd.read_csv('../Metadata.csv')
tkmers = pd.read_csv('../Top_kmers.csv', index_col = 0)

species = tkmers.columns.values[1:len(tkmers.columns.values)]
meta_ref = metadata[metadata['Number of Contigs'] == 2]
meta_ref = meta_ref.drop_duplicates(subset = 'Species', keep='first')

def refrence_location(species):
	ref_row = meta_ref[meta_ref['Species']==species]
	ref_file = ref_row['Sample'].to_list()
	ref_file_location = '../Approved_Sequences/'+ref_file[0]+'.fna'
	return(ref_file_location)

def refrence_file(species):
	location = refrence_location(species)
	sequence = []
	with open(location) as handle:
		for title, seq in SimpleFastaParser(handle):
			sequence.append(seq)
	chromosome1 = sequence[0]
	chromosome2 = sequence[1]
	return(chromosome1, chromosome2)

def needles(species):
	kmer_df = tkmers[species].dropna(axis='index', how='all')
	kmers = kmer_df.index.values
	aho = ahocorasick.Automaton()
	for i in range(0,len(kmers)):
		aho.add_word(kmers[i], (i, kmers[i]))
	aho.make_automaton()
	return(aho)

def find_needles(species):
	aho = needles(species)
	print(aho)
	#print(aho.get_stats())
	#print(aho.get('ACGTAGCCCGATCGCTATTGCGCTTTTCATG'))
	chromosome1, chromosome2 = refrence_file(species)
	it = aho.iter(chromosome1)
	print(it)
	for index, value in it:
		print(index, value) 
	for item in aho.iter(chromosome1):
		print("EXCUSE MEEEE")


find_needles('Brucella sp')