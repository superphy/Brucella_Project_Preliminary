def file_locations():
	directories = (os.listdir("refseq/bacteria")) # list of all the folders in the bacteria folder
	internal_directories = [] # list of the files in each folder
	internal_directories_locations = [] 
	fna_files = [] # list of the fna files
	fna_file_location = [] #the location of the fna files

	for file in directories:
		location_1 = "refseq/bacteria/"+file
		internal_directories.append(os.listdir(location_1))
		internal_directories_locations.append(location_1)

	while i< len(internal_directories):
		first_entry = internal_directories[i][0]
		file_location =""
		if first_entry == "MD5SUMS":
			fna_files.append(internal_directories[i][1])
			file_location = internal_directories_locations[i]+"/"+internal_directories[i][1]
			fna_file_location.append(file_location)
		else:
			fna_files.append(first_entry)
			file_location = internal_directories_locations[i]+"/"+first_entry
			fna_file_location.append(file_location)
		i=i+1
	return fna_file_location