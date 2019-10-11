import os 
import pandas as pd 
#IF YOU EDIT THIS FILE SNAKEMAKE WILL TRY TO RE-RUN KSNP

original_files = os.listdir("Approved_Sequences/")
input_list = []
cwd = os.getcwd()

#formatting
for file in original_files: 
	entry = cwd+"/Approved_Sequences/"+file+"	"+file[4:13]+'\n'
	input_list.append(entry)

#writing output file
with open("kSNP3_input.txt", "w") as output:
	for row in input_list:
		output.write(str(row))
