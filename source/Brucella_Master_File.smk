
rule all:
	input: "quast_output/quast_complete"

rule ncbi_data_retrieval:
	output:
		"refseq/"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'

rule metadata_creation:
	input:
		flag = "refseq/",
		un = "source/unzip.py",
		met = "source/metadata.py"
	output:
		"Metadata.csv"
	run:
		shell("python {input.un}")
		shell("python {input.met}")

rule quast:
	input:
		md = "Metadata.csv", 
		qu = "source/quast.py"
	output:
		"quast_output/quast_complete"	
	run:
		shell("python {input.qu}")