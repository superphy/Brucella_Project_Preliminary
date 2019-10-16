import os 
import pandas as pd 
i=0
metadata = pd.read_csv('Metadata-v1.csv')
samples = metadata['Sample'].tolist()

while i< len(samples):
	sample = samples[i]
	location = 'refseq/bacteria/'+sample[0:15]+'/'+sample+'.fna'
	sysin = "python quast-5.0.2/quast.py --silent --fast -o quast_files/quast_"+sample+' '+location
	os.system(sysin)
	if i == len(samples)-1:
		os.system("mkdir quast_files/quast_comp")
		print("quast complete")
	i+=1


