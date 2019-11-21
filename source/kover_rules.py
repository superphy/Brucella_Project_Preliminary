import pandas as pd
from Bio import SeqIO
import os.path
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

strain_occ = pd.read_csv('Strain_Occurances.csv', index_col = 0)
strains = strain_occ.index.values.tolist()
strains.remove('Brucella neotomae')
strains.remove('Brucella sp')
strains = [strain.replace(' ', '_') for strain in strains]

def rule_reader (rule): #Reads the fasta file of the given rule (if it exists) and returns a list of all the sequences and what the rule is
	rules = [str() for c in 'c' * 10000] #Create a list of empty 10000 strings (the number of possible rules)
	header = ''
	i = 0 
	if os.path.isfile(rule): #If the rule exists exists
		for record in SeqIO.parse(rule, 'fasta'):
			seq = record.seq
			seq = seq._get_seq_str_and_check_alphabet(seq)
			rules[i]=seq #Place the sequence into the appropriate location of the list of empty strings
			header = record.id 
			i+=1 
	return(header, rules)

def rules_merger(strain): #Calls Rule reader for all the possible rules
	rule1 = 'Kover_Data/'+strain+'/model_rule_1_equiv.fasta'
	rule2 = 'Kover_Data/'+strain+'/model_rule_2_equiv.fasta'
	header1, rules1 = rule_reader(rule1)
	header2, rules2 = rule_reader(rule2)
	return(strain, header1, rules1, header2, rules2)

df = pd.DataFrame()

with ProcessPoolExecutor(len(strains)) as ppe:
	for strain, header1, rules1, header2, rules2 in ppe.map(rules_merger, strains):
		df[strain+'-'+header1]=rules1
		if header2 != '':
			df[strain+'-'+header2]=rules2

df.to_csv('Kover_Rules.csv', index = False)