
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
		"Metadata.csv"
	output:
		"quast_output/quast_complete"	
	run:
		import pandas as pd 
		i=0
		metadata = pd.read_csv('Metadata.csv', sep='\t')
		print(metadata)
		fna_file_locations = metadata['File Location'].tolist()
		fna_file_names = metadata['File Name'].tolist()

		os.system("mkdir quast_output")

		while i< len(fna_file_locations):
			sys_in = "python quast-5.0.2/quast.py -o quast_output/quast_"+fna_file_names[i]+" "+fna_file_locations[i]
			os.system(sys_in)
			if i == len(fna_file_locations)-1:
				os.system("mkdir quast_output/quast_complete")
				print("quast complete")
			i = i+1
