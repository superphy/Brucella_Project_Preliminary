rule ncbi_data_retreval:
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria '
