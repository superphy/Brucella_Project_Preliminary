import os 
import pandas as pd 
import matplotlib.pyplot as plt

contigs_list = []
contigs_under_1000_list =[]
j=0

metadata = pd.read_csv('Metadata.csv')

while j<len(metadata):
	sample = metadata['Sample'].iat[j]
	report_location = "quast_files/quast_"+sample+"/report.tsv"
	report = pd.read_csv (report_location, sep = '\t')
	total_contigs = report.iat[12,1]
	contigs_list.append(total_contigs)
	contigs_under_1000 = total_contigs - report.iat[1,1]
	contigs_under_1000_list.append(contigs_under_1000)
	j+=1

metadata['Number of Contigs']= contigs_list
metadata['Number of Contigs <1000 bp'] = contigs_under_1000_list

metadata.to_csv('Metadata.csv', index = False )

scatterplot = plt.scatter(contigs_list, contigs_under_1000_list, s=10, c='k')
plt.xlabel ("Number of Contigs")
plt.ylabel("Number of Contigs <1000bp")
 
plt.title('Data Quality Visualization')
plt.savefig("Data Quality Visualization.png")
