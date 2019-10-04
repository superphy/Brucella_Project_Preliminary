import os 
import pandas as pd

files = os.listdir("refseq/bacteria")
fna_files_rm = ""
fna_files = []
n=5

# Deletes all the MD5SUMS files
for file in files:
	sysin = "rm refseq/bacteria/"+file+"/MD5SUMS"
	os.system(sysin)

# gathers the name of all the fna files to be run in refseq_masher
for file in files:
	sysin = "refseq/bacteria/"+file
	fna_files.append(os.listdir(sysin)[0])
	fna_files_rm= fna_files_rm+" "+sysin+"/"+os.listdir(sysin)[0]

#runs refseq_masher 
sysin = "refseq_masher matches "+fna_files_rm+" -o rmash_output.csv /home/ashlynn/Desktop/Fall_2019/random --output-type csv -n "+str(n)
os.system(sysin)