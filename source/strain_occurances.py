import os 
import pandas as pd
import re
sp_dict = {}
out_list = []
species = []
j=0
n=5

metadata = pd.read_csv("Metadata.csv")
strains = ['Brucella abortus','Brucella canis','Brucella ceti', 'Brucella melitensis', 'Brucella neotomae', 'Brucella ovis', 'Brucella pinnipedialis', 'Brucella sp', 'Brucella suis']

for key in strains:
	sp_dict[key] = 0

while j< len(metadata):
	#incrementing the correct counter for every occurance of a species name within a group 
	tax_sp = metadata['Species'].iat[j]
	match= re.search('Brucella sp.',tax_sp)
	if match != None: 
		tax_sp_ver = 'Brucella sp'
	if match == None:
		tax_sp_ver = tax_sp
	sp_dict[tax_sp_ver]+=1
	j+=1

df = pd.DataFrame.from_dict(sp_dict, orient = 'index', columns = ['Number of Occurances'])
df.to_csv('Strain_Occurances.csv')
