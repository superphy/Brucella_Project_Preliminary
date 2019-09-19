import os 
import pandas as pd 

directories = (os.listdir("quast_output/"))
quast_summary =[]
i=0

while i < len(directories):
	if directories[i] != "quast_complete":
		quast_file_name = directories[i]
		report_location = "quast_output/"+quast_file_name+"/report.tsv"
		report = pd.read_csv(report_location, sep='\t')
		data = report.iloc[:,1].to_list()
		contigs = data[12]
		contigs_less_1000 = contigs - data[1]
		row = [quast_file_name, contigs, contigs_less_1000]
		quast_summary.append(row)
	i=i+1

df =pd.DataFrame.from_records(quast_summary, columns =['File Name', 'Number of Contigs', 'Number of Contigs < 1000bp'])
df.to_csv('Quast_Summary.csv', index=False)

