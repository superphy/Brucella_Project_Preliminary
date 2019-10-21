fasta_files = "../Approved_Sequences/"
kmer_size = 31 

ids, = glob_wildcards(fasta_files+"{id}.fna")

rule all: 
	input:
		expand("mer_counts/{id}.fa", id = ids)

rule count:
  input:
    fasta_files+"{id}.fna"
  output:
    temp("mer_counts/{id}.jf")
  threads:
    2
  shell:
    "jellyfish count -C -m {kmer_size} -s 100M -t {threads} {input} -o {output} --out-counter-len 1"

rule dump:
    input:
        "mer_counts/{id}.jf"
    output:
        mc = "mer_counts/{id}.fa", 
    shell:
        "jellyfish dump {input} > {output.mc}"