import os 
import pandas as pd 
i=0
metadata = pd.read_csv('Metadata.csv', sep='\t')
fna_file_locations = metadata['File Location'].tolist()
fna_file_names = metadata['File Name'].tolist()

os.system("mkdir quast_output")

while i< len(fna_file_locations):
	sys_in = "python quast-5.0.2/quast.py -o quast_output/quast_"+fna_file_names[i]+" "+fna_file_locations[i]
	os.system(sys_in)
	if i == len(fna_file_locations)-1:
		os.system("mkdir quast_output/quast_complete")
		print("quast complete")
	i = i+1
