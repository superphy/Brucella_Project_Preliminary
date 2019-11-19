import pandas as pd
from Bio import SeqIO
import os.path
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

strain_occ = pd.read_csv('../Strain_Occurances.csv', index_col = 0)
strains = strain_occ.index.values.tolist()
strains.remove('Brucella neotomae')
strains.remove('Brucella sp')
strains = [strain.replace(' ', '_') for strain in strains]

def rules_parser(strain):
	rule1 = '../Kover_Data/'+strain+'/model_rule_1_equiv.fasta'
	rule2 = '../Kover_Data/'+strain+'/model_rule_2_equiv.fasta'
	rules1, rules2 = [str() for c in 'c' * 10000], [str() for c in 'c' * 10000]
	header1, header2 = '', ''
	i1, i2 = 0,0
	if os.path.isfile(rule1):
		for record in SeqIO.parse(rule1, 'fasta'):
			seq = record.seq
			seq = seq._get_seq_str_and_check_alphabet(seq)
			rules1[i1]=seq
			header1 = record.id
			i1+=1 
	if os.path.isfile(rule2):
		for record in SeqIO.parse(rule2, 'fasta'):
			seq = record.seq
			seq = seq._get_seq_str_and_check_alphabet(seq)
			rules2[i2]=seq
			header2 = record.id
			i2+=1
	return(strain, header1, rules1, header2, rules2)

df = pd.DataFrame()

with ProcessPoolExecutor(len(strains)) as ppe:
	for strain, header1, rules1, header2, rules2 in ppe.map(rules_parser, strains):
		df[strain+'-'+header1]=rules1
		if header2 != '':
			df[strain+'-'+header2]=rules2

df.to_csv('../Kover_Rules.csv', index = False)