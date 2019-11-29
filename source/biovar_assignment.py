import pandas as pd 

bpm = pd.read_csv('../Blast_Primer_Matching.csv')
files = list(bpm['File ID'].unique())

biovar_1 = [773, 774, 775, 424, 425, 426, 196, 197, 198]
biovar_2 = [773, 774, 775, 550, 551, 552, 277, 278, 279]
biovar_3 = [773, 774, 775, 298, 299, 300, 196, 197, 198]
biovar_4 = [773, 774, 775, 613, 614, 615, 196, 197, 198]
biovar_5 = [773, 613, 277, 196]

def biovar(file):
	file_df = bpm[bpm['File ID']==file].fillna(0)
	p1 = file_df['Primer One (BMEI1426-BMEI1427)'].to_list()
	p2 = file_df['Primer Two (BR1080f-BR1080r)'].to_list()
	p3 = file_df['Primer Three (BMEI1688-BMEI1687)'].to_list()
	p4 = file_df['Primer Four (BMEI0205f-BMEI0205r)'].to_list()

	amplicons = p1+p2+p3+p4
	amplicons = [int(amp) for amp in amplicons if amp != 0]

	if len(amplicons)<3:
		return('No Matched Biovar')
	else:
		bv1 = all(amps in biovar_1 for amps in amplicons)
		bv2 = all(amps in biovar_2 for amps in amplicons)
		bv3 = all(amps in biovar_3 for amps in amplicons)
		bv4 = all(amps in biovar_4 for amps in amplicons)
		bv5 = all(amps in amplicons for amps in biovar_5)

		if bv1:
			return('Biovar 1')
		if bv2:
			return('Biovar 2')
		if bv3:
			return('Biovar 3')
		if bv5:
			return('Biovar 5')
		if bv4:
			return('Biovar 4')
		else:
			return('No Matched Biovar')

df = pd.DataFrame()
for file in files:
	df[file] = [biovar(file)]

df = df.T
df.to_csv('../Suis_Biovar.csv')