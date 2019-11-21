ids1, = glob_wildcards("../Kover_Data/Kover_{id1}_Metadata.tsv")

rule all:
	input:
		expand("../Kover_Data/{id1}/results.json", id1 = ids1)

rule kover_create:
	input:
		"../Kover_Data/Kover_{id1}_Metadata.tsv"
	output:
		"../Kover_Data/{id1}/{id1}_dataset.kover"
	run:
		#export PATH=/home/ashlynn/Desktop/Fall_2019/Brucella/kover/bin/:$PATH
		shell('kover dataset create from-tsv --genomic-data ../Kover_Data/Kmer_Matrix.tsv --phenotype-description "{wildcards.id1} vs non {wildcards.id1}" --phenotype-metadata {input} --output {output} --progress ')


rule kover_split:
	input:
		"../Kover_Data/{id1}/{id1}_dataset.kover"
		#expand("../Kover_Data/{id1}_dataset.kover", id1 = ids1)
	output:
		temp('../Kover_Data/{id1}/{id1}_flag.txt')
	run:
		shell('kover dataset split --dataset {input} --id {wildcards.id1} --train-size 0.666 --folds 5 --random-seed 72 --progress'),
		shell('touch ../Kover_Data/{wildcards.id1}/{wildcards.id1}_flag.txt')

rule kover_learn:
	input:
		ds = "../Kover_Data/{id1}/{id1}_dataset.kover", 
		split_flag = '../Kover_Data/{id1}/{id1}_flag.txt'
	output:
		'../Kover_Data/{id1}/results.json'
	run:
		shell('kover learn scm --dataset {input.ds} --split {wildcards.id1} --hp-choice bound --max-rules 5  --progress --output-dir ../Kover_Data/{wildcards.id1}/')

'''
TEST

Create: 
	kover dataset create from-tsv --genomic-data Kmer_Matrix.tsv --phenotype-description "test" --phenotype-metadata Kover_Brucella_abortus_Metadata.tsv --output test_dataset.kover --progress 
Split:
	kover dataset split --dataset test_dataset.kover --id test --train-size 0.666 --folds 5 --random-seed 72 --progress
Learn:
	kover learn scm --dataset test_dataset.kover --split test --progress 
	
* Kover tests folder with matrix, metadata and dataset all together runs correctly - also works from one directory above & with specified output dir 
* WORKED - kover learn scm --dataset ../Kover_Data/Brucella_abortus/Brucella_abortus_dataset.kover --split Brucella_abortus --progress --output-dir ../Kover_Data/Brucella_abortus/
	run directly 
* FAILED - kover learn scm --dataset ../Kover_Data/Brucella_pinnipedialis/Brucella_pinnipedialis_dataset.kover --split Brucella_pinnipedialis  --progress --output-dir ../Kover_Data/Brucella_pinnipedialis/ 
	run through snakemake & directly
	'''