from multiprocessing import cpu_count
kmer_size = 31
cpu = cpu_count() - 6


rule all:
	input:
		'j_all_brucella.fna'

#Aquiring the data from the ncbi database
rule ncbi_data_retrieval:
	output:
		"refseq/"
	shell:
		'ncbi-genome-download --parallel 55 --genus "brucella" bacteria --format fasta'

#unzipping the downloaded files
rule unzip: 
	input:
		ref = 'refseq/',
		un = 'source/unzip.py'
	output:
		'refseq/dunzip'
	run:
		shell('python {input.un}')

#Running refseq_masher on the NCBI data
rule rmash:
	conda:
		'../envs/mash.yaml'
	input:
		rm = 'source/r_mash.py',
		ref = "refseq/"
	output: 
		'rmash_output.csv'
	shell:
		'python {input.rm}'

#Compressing the output from rmash
rule rmash_filtering:
	input:
		rmo = 'rmash_output.csv',
		rmf = 'source/r_mash_filtering.py'
	output: 
		'rmash_filtered_output.csv'
	run: 
		shell('python {input.rmf}')

# Comparing the species listed in the fasta file and the results of rmash
rule species_val:
	input: 
		un = 'refseq/dunzip', 
		sv = 'source/species_validation.py'
	output: 
		'Metadata-v1.csv'
	run:
		shell('python {input.sv}')
'''
#Running quast on all of the datasets
rule quast:
	input:
		md = "Metadata-v1.csv", 
		qu = "source/quast.py"
	output:
		"quast_files/quast_comp"
	run:
		shell("python {input.qu}")

#Summarizing the results of the analysis completed by quast 
rule quast_summary:
	input:
		flag = "quast_files/quast_comp",
		qs = "source/quast_summary.py"
	output:
		"Metadata-v2.csv"
	run:
		shell("python {input.qs}")
'''
#Removing any contigs with under 500bp
rule contig_len_filtering:
	input:
		md = "Metadata-v2.csv",
		clf = "source/contig_len_filtering.py"
	output:
		md = "Metadata.csv",
		aps = 'Approved_Sequences/'
	run:
		shell("python {input.clf}")		

#Generates a csv with the breakdown of samples per species in the dataset
rule strain_occurances: 
	input: 
		md = "Metadata.csv",
		so = "source/strain_occurances.py"
	output:
		"Strain_Occurances.csv"
	run:
		shell("python {input.so}")


#Generating the input file for kSNP3
rule ksnp_input:
	input:
		#aps = 'Approved_Sequences/',
		ki = "source/ksnp_input.py"
	output: 
		"kSNP3_input.txt"	
	run:
		shell("python {input.ki}")

#Running kSNP3
rule ksnp:
	input:
		"kSNP3_input.txt"
	output:
		"kSNP3_Output/Logfile.txt"
	run:
		#shell("PATH=$PATH:~/Desktop/Fall_2019/Brucella/kSNP3.1_Linux_package/kSNP3")
		shell("kSNP3 -in kSNP3_input.txt -outdir kSNP3_Output/ -k 31 -CPU {cpu} | tee kSNP3_Output/Logfile.txt")

#Building the tree from the Newick file
rule tree:
	input:
		md = 'Metadata.csv',
		ko = "kSNP3_Output/Logfile.txt",
		tm = "source/tree_manipulation.py"
	output:
		"tree.pdf"
	run:
		shell("python {input.tm}")

# Generates a fasta file that is the sum of all fasta files in the dataset
rule all_brucella:
	input:
		'Approved_Sequences/'
	output:
		'all_brucella.fna'
	run:
		shell('cat {input}*.fna > {output}')

rule all_jellyfish_count: 
	input: 
		'all_brucella.fna'
	output:
		temp('j_all_brucella.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t {cpu} {input} -o {output} -L 190 -U 501 --out-counter-len 1')

rule all_jellyfish_dump: 
	input: 
		'j_all_brucella.jf'
	output:
		dump = 'j_all_brucella.fna'
	run:
		shell('jellyfish dump {input} > {output}')