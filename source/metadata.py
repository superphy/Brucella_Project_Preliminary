import os 
from Bio import SeqIO
import re 
import csv
import pandas as pd 

directories = (os.listdir("refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories = [] # list of the files in each folder
internal_directories_locations = [] 
fna_files = [] # list of the fna files
fna_file_location = [] #the location of the fna files
sequence_data = [] #list of records

i=0 
j=0

for file in directories:
	location_1 = "refseq/bacteria/"+file
	internal_directories.append(os.listdir(location_1))
	internal_directories_locations.append(location_1)

while i< len(internal_directories):
	first_entry = internal_directories[i][0]
	file_location =""
	if first_entry == "MD5SUMS":
		fna_files.append(internal_directories[i][1])
		file_location = internal_directories_locations[i]+"/"+internal_directories[i][1]
		fna_file_location.append(file_location)
	else:
		fna_files.append(first_entry)
		file_location = internal_directories_locations[i]+"/"+first_entry
		fna_file_location.append(file_location)
	i=i+1

while j < len(fna_file_location):
	for seq_record in SeqIO.parse(fna_file_location[j], "fasta"):
		record_id =[]
		match_a = re.search('N\w\w\w\w\w\w',seq_record.id)
		if match_a != None:
			id_short = match_a.group(0)
		if match_a == None:
			print(seq_record.id)
			id_short = "no id found"
		record_id.append(id_short)
		
		description = seq_record.description
		match_b = re.search('Brucella\s\w+',description)
		if match_b != None:
			strain = match_b.group(0)
		if match_b == None:
			strain = "no strain found"	
	record = [fna_file_location[j], fna_files[j], record_id[0], strain]
	sequence_data.append(record)
	j=j+1

df =pd.DataFrame.from_records(sequence_data, columns =['File Location','File Name','Record ID','Strain'])
df.to_csv('Metadata.csv', index = False)

print ("Done creating Metadata.csv")