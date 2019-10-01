import pandas as pd 
import os
from Bio import SeqIO
i=0

quast_summary = pd.read_csv('Quast_Summary.csv')
metadata = pd.read_csv('Metadata.csv')
files_w_contigs_under_1000bp =[]
files_wo_contigs_under_1000bp =[]
fna_names = []
fna_file_locations = []
file_location_list = []
file_location_list_b = []

while i< len(quast_summary):
	quast_row = quast_summary.iloc[i]
	quast_name = quast_row.iloc[1]
	fna_location = quast_row.iloc[0]
	fna_name = quast_name[6:len(quast_name)]
	file_info = [fna_location,fna_name,quast_row.iloc[1],quast_row.iloc[2]]
	if quast_row.iloc[3] > 0:
		files_w_contigs_under_1000bp.append(file_info)	
	else: 
		files_wo_contigs_under_1000bp.append(file_info)	
	i=i+1

#os.system("mkdir Approved_Sequences")

for record in files_wo_contigs_under_1000bp:
	cp = "cp "+record[0]+" Approved_Sequences/"
	fll = "Approved_Sequences/"+record[1]
	file_location_list.append(fll)
	os.system(cp)

for record in files_w_contigs_under_1000bp:
	len_500_plus =[]
	input_location = record[0]
	for seq_record in SeqIO.parse(input_location, "fasta"):
		if len(seq_record.seq)>500:
			len_500_plus.append(seq_record)
	output_location = "Approved_Sequences/"+record[1]
	file_location_list.append(output_location)		
	SeqIO.write(len_500_plus,output_location, "fasta")	


print ("Done contig length filtering")	

