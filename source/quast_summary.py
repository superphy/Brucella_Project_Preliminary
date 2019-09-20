import os 
import pandas as pd 
import matplotlib.pyplot as plt

directories = (os.listdir("quast_output/"))
contigs_list = []
contigs_less_1000_list =[]
quast_summary =[]
i=0

while i < len(directories):
	if directories[i] != "quast_complete":
		quast_file_name = directories[i]
		report_location = "quast_output/"+quast_file_name+"/report.tsv"
		report = pd.read_csv(report_location, sep='\t')
		data = report.iloc[:,1].to_list()
		contigs = data[12]
		contigs_list.append(contigs)
		contigs_less_1000 = contigs - data[1]
		contigs_less_1000_list.append(contigs_less_1000)
		og_file_location = "refseq/bacteria/"+quast_file_name[6:21]+"/"+quast_file_name[6:len(quast_file_name)]
		row = [og_file_location, quast_file_name, contigs, contigs_less_1000]
		quast_summary.append(row)
	i=i+1

df =pd.DataFrame.from_records(quast_summary, columns =['Fna File Location','File Name', 'Number of Contigs', 'Number of Contigs < 1000bp'])
df.to_csv('Quast_Summary.csv', index=False)

scatterplot = plt.scatter(contigs_list, contigs_less_1000_list, s=10, c='k')
plt.xlabel ("Number of Contigs")
plt.ylabel("Number of Contigs <1000bp")
 
plt.title('Data Quality Visualization')
plt.savefig("Data Quality Visualization.png")

