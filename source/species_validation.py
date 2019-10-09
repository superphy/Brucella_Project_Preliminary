import os 
import pandas as pd
from Bio import SeqIO
import re 

i=0
j=0

rmash = pd.read_csv("../rmash_filtered_output.csv")
new_col = []

while i<len(rmash):
	sample = rmash['Sample'].iat[i]
	location = "../refseq/bacteria/"+sample[0:15]+"/"+sample+".fna"
	for seq_record in SeqIO.parse(location, "fasta"):
		description = seq_record.description
		match= re.search('Brucella\s\w+',description)
		if match!= None:
			fna_strain = match.group(0)
		if match== None:
			fna_strain = "no strain found"	
	new_col.append(fna_strain)
	i+=1
rmash.insert(loc = 2, column='Fna Species', value = new_col)

diff_results = 0

while j<len(rmash):
	if rmash['Fna Species'].iat[j] == 'Brucella sp':
		rmash['Fna Species'].iat[j] = rmash['Rmash Species'].iat[j]

	if rmash['Fna Species'].iat[j] != rmash['Rmash Species'].iat[j]:
		diff_results +=1
	j+=1 

print(diff_results)