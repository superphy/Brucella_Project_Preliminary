rule ncbi_data_retrieval:
	conda: 
		"envs/ncbi_file_retrieval.yaml"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'
