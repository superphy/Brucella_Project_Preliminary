
rule all:
	input: "Approved_Sequences/"

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

rule quast_summary:
	input:
		flag = "quast_output/quast_complete",
		qs = "source/quast_summary.py"
	output:
		"Quast_Summary.csv"
	run:
		shell("python {input.qs}")

rule contig_len_filtering:
	input:
		qs = "Quast_Summary.csv", 
		md = "Metadata.csv",
		clf = "source/contig_len_filtering.py"
	output:
		"Approved_Sequences/"
	run:
		shell("python {input.clf}")		
