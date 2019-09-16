import os 
import Bio
from Bio import SeqIO

directories = (os.listdir("refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories_two = [] # list of the files in each folder
internal_directories_one = [] # list of the fna files
fna_file_paths =[] # list of the paths to all the fna files
i=0 
j=0

for file in directories:
	location_1 = "refseq/bacteria/"+file
	internal_directories_two.append(os.listdir(location_1))

while i< len(internal_directories_two):
	internal_directories_one.append(internal_directories_two[i][0]) 
	i=i+1

while j <len(internal_directories_two): 
	location_2 = "refseq/bacteria/"+directories[j]+"/"+internal_directories_one[j]
	fna_file_paths.append(location_2)
	j=j+1
	for seq_record in SeqIO.parse(location_2, "fasta"):
		print(seq_record.id)
		print(repr(seq_record.seq))
		print(len(seq_record))
		print("next")


