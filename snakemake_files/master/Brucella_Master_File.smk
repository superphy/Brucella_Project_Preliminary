rule ncbi_data_retrieval:
	conda: 
		"../envs/ncbi_file_retrieval.yaml"
	output:
		"../.."
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'

rule unzip_files: 
	input:
		"../../python_files/unzip.py" 
	shell:
		"python unzip.py"
