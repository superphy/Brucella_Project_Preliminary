import os 

directories = (os.listdir("refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories_two = [] # list of the files in each folder
internal_directories_one = [] # list of the fna files
i =0 

for file in directories:
	location = "refseq/bacteria/"+file
	internal_directories_two.append(os.listdir(location))

while i< len(internal_directories_two):
	internal_directories_one.append(internal_directories_two[i][0]) 
	i=i+1

print(internal_directories_one)
print(len(internal_directories_one))
