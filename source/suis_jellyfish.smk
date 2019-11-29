import pandas as pd 
metadata = pd.read_csv('../Metadata.csv')
suis_meta = metadata[metadata['Species']=='Brucella suis']
samples = list(suis_meta['Sample'])
samples = [sample+'.fna' for sample in samples]

#Performs the jellyfish count opperation on all_brucella.fna
rule all_jellyfish_count: 
	input: 
		'../blast/suis_master.fna'
	output:
		temp('j_all_suis.jf')
	run:
		shell('jellyfish count -C -m {kmer_size} -s 100M -t {cpu} {input} -o {output} -L 190 -U 501 --out-counter-len 1')

#Performs the jellyfish dump opperation on all_brucella.fna
rule all_jellyfish_dump: 
	input: 
		'j_all_suis.jf'
	output:
		'j_all_suis.fna'
	run:
		shell('jellyfish dump {input} > {output}')

#Performs the jellyfish count opperation on everything in the dataset
rule count:
	input:
		"../Approved_Sequences/{samples}.fna"
	output:
		temp("../suis_jellyfish_output/{samples}.jf")
	shell:
		"jellyfish count -C -m {kmer_size} -s 100M -t 2 {input} -o {output} --out-counter-len 1"

#Performs the jellyfish dump opperation on everything in the dataset
rule dump:
	input:
		"../suis_jellyfish_output/{samples}.jf"
	output:
		"../suis_jellyfish_output/{samples}.fna", 
	shell:
		"jellyfish dump {input} > {output}"

#Makes count/dump actually work
rule jellyfish_complete:
	input:
		expand("../suis_jellyfish_output/{samples}.fna", samples = samples)
	output:
		'../flags/suis_jellyfish_flag.txt'
	run:
		shell('touch ../flags/jellyfish_flag.txt')
'''
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
'''