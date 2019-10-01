import os 
import pandas as pd 

files = os.listdir("Approved_Sequences/")
input_list = []
cwd = os.getcwd()


for file in files:
	file_location = cwd+"/"+file
	id_name = file[4:13]
	entry = file_location+"	"+id_name+'\n'
	input_list.append(entry)

with open("kSNP3_input.txt", "w") as output:
	for row in input_list:
		output.write(str(row))

