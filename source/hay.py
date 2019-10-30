#!/usr/bin/env python

# Example output:
# C:\Users\LaingC\dev\Brucella\brucella-env\Scripts\python.exe C:/Users/LaingC/dev/Brucella/source/hay.py ..\Approved_Sequences\GCF_000007125_1_ASM712v1_genomic.fna
# NC_003317.1 Brucella melitensis bv. 1 str. 16M chromosome I, complete sequence
# 30 TAATGATAAATGATAACACAGGATGGAATGG
# 2132 TCCTTATATTTAAAAAGAGTTCCTCATCCAT
# 10134 AATGCAATGGCCGCTGCCGACTCCGTTCTGG
# NC_003318.1 Brucella melitensis 16M chromosome II, complete sequence
#
# Process finished with exit code 0


import sys
import ahocorasick
from Bio import SeqIO

# get filename from commandline
genome_fasta = sys.argv[1]

a_kmers = ahocorasick.Automaton()

# add some kmers you want to look for -- in reality these would
# be from an external source, but these are for example
# be sure to account for the reverse complement eg. revcomp(kmer) in the actual code
kmers_from_file = ["TAATGATAAATGATAACACAGGATGGAATGG", "TCCTTATATTTAAAAAGAGTTCCTCATCCAT", "AATGCAATGGCCGCTGCCGACTCCGTTCTGG"]

# add the kmers to your automaton
# we want the kmer sequence as the value as well, so we can access it after a match
for k in kmers_from_file:
    a_kmers.add_word(k, k)

a_kmers.make_automaton()

# search the genome for the kmers
with open(genome_fasta) as gfile:
    for record in SeqIO.parse(gfile, "fasta"):
        #print out the fasta header
        print(record.description)

        # search each contig for the kmers
        # the results are the end index of the match (the last bp of the kmer) and whatever value
        # we stored -- which was the kmer sequence itself
        for (end_bp, kmer) in a_kmers.iter(str(record.seq)):
            print(end_bp, kmer)


