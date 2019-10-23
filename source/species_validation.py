import os 
import pandas as pd
from Bio import SeqIO
import re 

i=0
j=0

rmash = pd.read_csv("rmash_filtered_output.csv")
fna_sp_col = []
rec_id_col = []
final_sp_col = []

#gathering metadata information and fna species info
while i<len(rmash):
	sample = rmash['Sample'].iat[i]
	location = "refseq/bacteria/"+sample[0:15]+"/"+sample+".fna"
	for seq_record in SeqIO.parse(location, "fasta"):
		
		description = seq_record.description
		match= re.search('Brucella\s\w+',description)
		if match!= None:
			fna_strain = match.group(0)
		if match== None:
			fna_strain = "no strain found"	
	
		record_id =[]
		match_a = re.search('N\w\w\w\w\w\w',seq_record.id)
		if match_a != None:
			id_short = match_a.group(0)
		if match_a == None:
			id_short = "no id found"
	rec_id_col.append(id_short)
	fna_sp_col.append(fna_strain)
	i+=1

# adding a Species and Record ID column to the dataframe 
rmash.insert(loc = 2, column='Species', value = fna_sp_col)
rmash.insert(loc = 3, column= 'Record ID', value = rec_id_col)


# comapring the rmash results and the fna species
while j<len(rmash):
	sample = rmash['Sample'].iat[j]
	# if the fna species is sp, it takes on the species found by rmash
	if rmash['Species'].iat[j] == 'Brucella sp':
		if rmash['Max Distance'].iat[j] < 0.001:
			rmash['Species'].iat[j] = rmash['Rmash Species'].iat[j]
		else:
			rmash['Rmash Species'].iat[j] = rmash['Species'].iat[j]
	# if the fna and rmash species info does not match the file is deleted  
	if rmash['Species'].iat[j] != rmash['Rmash Species'].iat[j]:
		sysin = 'rm -rf refseq/bacteria/'+sample[0:15]
		os.system(sysin)
	j+=1 

#deleting any entries where the fna and rmash species info does not match
rmash = rmash[rmash['Species']==rmash['Rmash Species']]

# deleting the Rmash column and creating metadata.csv 
rmash.drop(['Rmash Species'], axis = 1, inplace = True)
rmash.to_csv('Metadata-v1.csv', index = False)