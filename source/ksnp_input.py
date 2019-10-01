import os 
import pandas as pd 
i=0
#IF YOU EDIT THIS FILE SNAKEMAKE WILL TRY TO RE-RUN KSNP

original_files = os.listdir("Approved_Sequences/")
input_list = []
cwd = os.getcwd()

for file in original_files:
	new_name = file[0:len(file)-4]
	newer_name = ""
	#replacing all the periods with underscores
	while i < len(new_name):
		if new_name[i] == ".":
			newer_name=new_name[0:i]+"_"+new_name[i+1:len(new_name)]
			break
		i=i+1
	#renaming the actual files	
	sysin = "mv Approved_Sequences/"+file+" Approved_Sequences/"+newer_name+".fna"
	os.system(sysin)
	entry = cwd+"/"+newer_name+".fna"+"	"+newer_name[4:13]+'\n'
	input_list.append(entry)

#writing output file
with open("kSNP3_input.txt", "w") as output:
	for row in input_list:
		output.write(str(row))
