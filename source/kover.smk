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
		'(kover learn scm --dataset ../Kover_Data/{wildcards.id1}/{wildcards.id1}_dataset.kover --split {wildcards.id1} --model-type conjunction disjunction --p 0.1 1.0 10.0 --max-rules 5 --hp-choice cv --n-cpu 2 --progress --output-dir ../Kover_Data/{wildcards.id1}/ )'
