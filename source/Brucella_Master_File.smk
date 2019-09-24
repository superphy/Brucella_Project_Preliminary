
rule all:
	input: "kSNP3_input.txt"

#Aquiring the data from the ncbi database
rule ncbi_data_retrieval:
	output:
		"refseq/"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'

#Creating a metadata file with the folder names and species
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

#Running quast on all of the datasets
rule quast:
	input:
		md = "Metadata.csv", 
		qu = "source/quast.py"
	output:
		"quast_output/quast_complete"
	run:
		shell("python {input.qu}")

#Summarizing the results of the analysis completed by quast 
rule quast_summary:
	input:
		flag = "quast_output/quast_complete",
		qs = "source/quast_summary.py"
	output:
		"Quast_Summary.csv"
	run:
		shell("python {input.qs}")

#Removing any contigs with under 500bp
rule contig_len_filtering:
	input:
		qs = "Quast_Summary.csv", 
		md = "Metadata.csv",
		clf = "source/contig_len_filtering.py"
	output:
		"Approved_Sequences/"
	run:
		shell("python {input.clf}")		

#Generating the input file for kSNP3
rule ksnp_input:
	input:
		md = "Metadata.csv", 
		ki = "source/ksnp_input.py"
	output: 
		"kSNP3_input.txt"	
	run:
		shell("python {input.ki}")	