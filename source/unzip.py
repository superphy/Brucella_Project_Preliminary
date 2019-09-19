import os 
import gzip

directories = (os.listdir("refseq/bacteria")) # list of all the folders in the bacteria folder
internal_directories_two = [] # list of the files in each folder
internal_directories_one = [] # list of the fna.gz files
i=0 
j=0

for file in directories:
	location_1 = "refseq/bacteria/"+file
	internal_directories_two.append(os.listdir(location_1))

while i< len(internal_directories_two):
	first_entry = internal_directories_two[i][0]
	if first_entry == "MD5SUMS":
		internal_directories_one.append(internal_directories_two[i][1])
	else:
		internal_directories_one.append(first_entry)
	i=i+1

while j <len(internal_directories_two): 
	location_2 = "refseq/bacteria/"+directories[j]+"/"+internal_directories_one[j]
	argument = "gunzip " + location_2
	os.system(argument)
	j=j+1

print("done unzipping")