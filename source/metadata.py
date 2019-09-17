import os 
from Bio import SeqIO
import re 
import csv

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
		description = seq_record.description
		match = re.search('Brucella\s\w+',description)
		if match != None:
			strain = match.group(0)
		if match == None:
			strain = "no strain found"	
		record = [fna_files[j], seq_record.id, strain]
		sequence_data.append(record)
	j=j+1

with open ('Metadata.csv', 'w') as output:
	output_writer = csv.writer(output, delimiter = '	')
	output.write('File Name, Record ID, Strain\n')
	for record in sequence_data:
		output_writer.writerow(record)