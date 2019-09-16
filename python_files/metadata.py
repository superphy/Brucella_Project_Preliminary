import os 
import gzip
import Bio 
from Bio import SeqIO

directories = (os.listdir("../refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories = [] # list of the files in each folder
internal_directories_locations = [] 
fna_files = [] # list of the fna files
fna_file_location = [] #the location of the fna files
i=0 
j=0

for file in directories:
	location_1 = "../refseq/bacteria/"+file
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

#while j < len(fna_file_location):
#	for seq_record in SeqIO.parse(fna_file_location[j], "fasta"):
#		print(seq_record.id)
#		print(seq_record.description)
		#print(repr(seq_record.seq))
#		print(len(seq_record))

while j < 1:
	for seq_record in SeqIO.parse(fna_file_location[j], "fasta"):
		fna_file =fna_files[j]
		description = seq_record.description
		#print(repr(seq_record.seq))
		sequence_lentgth = len(seq_record)
	j=j+1