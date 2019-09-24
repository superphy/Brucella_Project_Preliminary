import os 
import pandas as pd 

metadata = pd.read_csv('Metadata.csv')
input_list = []
i=0

while i<len(metadata):
	entry = metadata.iloc[i]
	k_entry = "Approved_Sequences/"+entry[1]+"	"+entry[2]+"\n"
	input_list.append(k_entry)
	i=i+1

with open("kSNP3_input.txt", "w") as output:
	for row in input_list:
		output.write(str(row))
