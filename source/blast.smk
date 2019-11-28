rule all:
	input:
		'Blast_Primer_Matching.csv'

rule blast_master_fasta:
	input:
		aps = 'Approved_Sequences/',
		sm = 'source/suis_meta.py'
	output:
		'blast/suis_master.fna'
	run:
		shell('python {input.sm}')

rule blast_database:
	input:
		master = 'blast/suis_master.fna'
	output:
		blast_db = 'blast/suis_master.fna.ndb'
	run:
		shell('makeblastdb -in {input} -parse_seqids -blastdb_version 5 -title "Brucella Suis" -dbtype nucl')

rule blast_search:
	input:
		query_file = 'queery_file.fna',
		blast_db = 'blast/suis_master.fna.ndb'
	output:
		serch_out = 'blast/blast_search_output.tsv'
	run:
		shell('blastn -db blast/suis_master.fna -query {input.query_file} -dust no -word_size 7 -evalue 100 -outfmt 6 -out {output}')

rule blast_output:
	input:
		blast_search_output = 'blast/blast_search_output.tsv',
		bo = 'source/blast_output.py'
	output:
		'Blast_Primer_Matching.csv'
	run:
		shell('python {input.bo}')