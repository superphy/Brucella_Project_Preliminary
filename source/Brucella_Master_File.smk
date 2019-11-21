from multiprocessing import cpu_count
kmer_size = 31
cpu = cpu_count() - 6

rule all:
	input:
		'Kover_Rules.csv'

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

#Generates a fasta file that is the sum of all fasta files in the dataset
rule all_brucella:
	input:
		'Approved_Sequences/'
	output:
		'all_brucella.fna'
	run:
		shell('cat {input}*.fna > {output}')

#Performs the jellyfish count opperation on all_brucella.fna
rule all_jellyfish_count: 
	input: 
		'all_brucella.fna'
	output:
		temp('j_all_brucella.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t {cpu} {input} -o {output} -L 190 -U 501 --out-counter-len 1')

#Performs the jellyfish dump opperation on all_brucella.fna
rule all_jellyfish_dump: 
	input: 
		'j_all_brucella.jf'
	output:
		'j_all_brucella.fna'
	run:
		shell('jellyfish dump {input} > {output}')

ids, = glob_wildcards("Approved_Sequences/{id}.fna")

#Performs the jellyfish count opperation on everything in the dataset
rule count:
	input:
		"Approved_Sequences/{id}.fna"
	output:
		temp("mer_counts/{id}.jf")
	shell:
		"jellyfish count -C -m {kmer_size} -s 100M -t 2 {input} -o {output} --out-counter-len 1"

#Performs the jellyfish dump opperation on everything in the dataset
rule dump:
	input:
		"mer_counts/{id}.jf"
	output:
		"mer_counts/{id}.fna", 
	shell:
		"jellyfish dump {input} > {output}"

#Makes count/dump actually work
rule jellyfish_complete:
	input:
		expand("mer_counts/{id}.fna", id = ids)
	output:
		'flags/jellyfish_flag.txt'
	run:
		shell('touch flags/jellyfish_flag.txt')

#Creates kmer_counts.csv 
rule kmer_dataframe:
	input:
		jf = 'flags/jellyfish_flag.txt',
		c_df = 'source/count_to_df.py'
	output:
		fmc = 'filtered_mer_counts/',
		kc = 'kmer_counts.csv'
	run:
		shell('python {input.c_df}')

#Analyzes the occurances of a specific kmer in each species
rule ranks:
	input:
		ra = 'source/ranks.py', 
		md = 'Metadata.csv', 
		so = 'Strain_Occurances.csv',
		kc = 'kmer_counts.csv'
	output:
		'Ranks.csv'
	run:
		shell('python {input.ra}')

#Generates histograms to visualize the rank distribution of each species
rule rank_hist:
	input:
		ranks = 'Ranks.csv',
		so = 'Strain_Occurances.csv',
		tk = 'source/rank_histograms.py'
	output:
		'Rank_Histograms/'
	run:
		shell('python {input.tk}')

#Locates each kmer on a refrence chromosome 
rule kmer_locations:
	input:
		cm = 'source/chrom_mapping.py', 
		md = 'Metadata.csv',
		ranks = 'Ranks.csv',
		so = "Strain_Occurances.csv"
	output:
		klv = 'Kmer_Location_Visualizations/', 
		kl = 'Kmer_Locations.csv'
	run:
		shell('python {input.cm}')

#Creates a folder of kover data (the kmer matrix and metadata files for each species)
rule kover_data:
	input:
		og_meta ='Metadata.csv',
		strain_occ = 'Strain_Occurances.csv',
		kcounts = 'kmer_counts.csv',
		ki = 'source/kover_inputs.py'
	output:
		kif = 'Kover_Data/',
		kmatrix = 'Kover_Data/Kmer_Matrix.tsv'
	run:
		shell('python {input.ki}')

ids1, = glob_wildcards("Kover_Data/Kover_{id1}_Metadata.tsv")

#Creates kover datasets from Metadata and Kmer Martrix (generated by kover_data)
rule kover_create:
	conda:
		'../envs/Kover.yaml'
	input:
		md = "Kover_Data/Kover_{id1}_Metadata.tsv",
		km = "Kover_Data/Kmer_Matrix.tsv"
	output:
		"Kover_Data/{id1}/{id1}_dataset.kover"
	shell:
		#export PATH=/home/ashlynn/Desktop/Fall_2019/Brucella/kover/bin/:$PATH
		'kover dataset create from-tsv --genomic-data Kover_Data/Kmer_Matrix.tsv --phenotype-description "{wildcards.id1} vs non {wildcards.id1}" --phenotype-metadata {input.md} --output {output} --progress '

#Splits the datasets into training and testing sets
rule kover_split:
	conda:
		'../envs/Kover.yaml'
	input:
		"Kover_Data/{id1}/{id1}_dataset.kover"
	output:
		temp('Kover_Data/{id1}/{id1}_flag.txt')
	shell:
		'kover dataset split --dataset {input} --id {wildcards.id1} --train-size 0.666 --folds 5 --random-seed 72 --progress && touch Kover_Data/{wildcards.id1}/{wildcards.id1}_flag.txt'

#Executes cover learn on the datasets to obtain markers for each dataset
rule kover_learn:
	conda:
		'../envs/Kover.yaml'
	input:
		ds = "Kover_Data/{id1}/{id1}_dataset.kover", 
		split_flag = 'Kover_Data/{id1}/{id1}_flag.txt'
	output:
		'Kover_Data/{id1}/results.json'
	shell:
		'kover learn scm --dataset {input.ds} --split {wildcards.id1} --hp-choice bound --max-rules 5  --progress --output-dir Kover_Data/{wildcards.id1}/'

#Makes the kover rules actually run
rule kover_complete:
	input:
		expand("Kover_Data/{id1}/results.json", id1 = ids1)
	output:
		temp('flags/kover_data_flag.txt')
	run:
		shell('touch flags/kover_data_flag.txt')

#Summarizes the rules established by Kover
rule kover_summary:
	input:
		flag = 'flags/kover_data_flag.txt',
		kr = 'source/kover_rules.py'
	output:
		'Kover_Rules.csv'
	run:
		shell('python {input.kr}')
