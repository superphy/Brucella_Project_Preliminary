import pandas as pd 
import os
from Bio import SeqIO
i=0
rename_list =[]

metadata = pd.read_csv('Metadata-v2.csv')

while i< len(metadata):
	output_name = metadata['Sample'].iat[i]
	location = "refseq/bacteria/"+output_name[0:15]+"/"+output_name+".fna"
	# Complying with kSNP3 naming standards 
	j=0
	while j < len(output_name):
		if output_name[j] == '.':
			output_name = output_name[0:j]+"_"+output_name[j+1:len(output_name)]
		j+=1

	output_location = "Approved_Sequences/"+output_name+".fna"
	# Removing contigs under 500 bp 
	if metadata['Number of Contigs <1000 bp'].iat[i] > 0:
		len_500_plus = []
		for seq_record in SeqIO.parse(location, "fasta"):
			if len(seq_record.seq)>500:
				len_500_plus.append(seq_record)		
		SeqIO.write(len_500_plus,output_location, "fasta")	
	# Copying the files without any contigs <1000bp to the Approved Sequences folder
	else:
		cp = "cp "+location+" Approved_Sequences/"	
		os.system(cp)
		mv = "mv Approved_Sequences/"+metadata['Sample'].iat[i]+".fna "+output_location
		rename_list.append(mv)
		
	metadata['Sample'].iat[i] = output_name
	i+=1

for mv in rename_list:
	os.system(mv)

metadata.to_csv('Metadata.csv', index = False)


