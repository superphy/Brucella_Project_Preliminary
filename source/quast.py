import os 
import pandas as pd 
i=0
metadata = pd.read_csv('Metadata.csv')
fna_file_locations = metadata['File Location'].tolist()
fna_file_names = metadata['File Name'].tolist()
quast_file_names = []

os.system("mkdir quast_output")

while i< len(fna_file_locations):
	file_name = "quast_"+fna_file_names[i]
	quast_file_names.append([file_name])
	sys_in = "python quast-5.0.2/quast.py -o quast_output/quast_"+fna_file_names[i]+" "+fna_file_locations[i]
	os.system(sys_in)
	i = i+1

df =pd.DataFrame.from_records(quast_file_names, columns =['File Name'])
df.to_csv('Quast_Summary.csv', index=False)

