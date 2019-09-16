import os 
import gzip
import Bio 
from Bio import SeqIO

directories = (os.listdir("../refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories = [] # list of the files in each folder
fna_files = [] # list of the fna files
i=0 
j=0

for file in directories:
	location_1 = "../refseq/bacteria/"+file
	internal_directories.append(os.listdir(location_1))

while i< len(internal_directories):
	first_entry = internal_directories[i][0]
	if first_entry == "MD5SUMS":
		fna_files.append(internal_directories[i][1])
	else:
		fna_files.append(first_entry)
	i=i+1

for file in fna_files:
	for seq_record in SeqIO.parse(file, "fasta"):
		print(seq_record.id)
		print(repr(seq_record.seq))
		print(len(seq_record))