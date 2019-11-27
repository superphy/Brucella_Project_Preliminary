import skbio.io
import pandas as pd
import numpy as np 
pd.options.mode.use_inf_as_na = True

blast_output_path = '../blast/blast_search_output.tsv'
with open(blast_output_path) as fh:
	blast_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True)

#blast_df = blast_df.set_index('qseqid')
blast_df = blast_df[blast_df['pident']==100]
#blast_df.to_csv('../Blast_Output.csv')


sseqids = list(blast_df['sseqid'].unique())

def min_max(start_a, end_a, start_b, end_b):
	primer_min = min(start_a, start_b)
	if primer_min == start_a:
		amp_start = max(start_a, end_a) 
		amp_end = min(start_b, end_b)
	if primer_min == start_b:
		amp_start = max(start_b, end_b) 
		amp_end = min(start_a, end_a)
	return(amp_start, amp_end)
'''
def primer_one(df, sseqid):
	BMEI1427 = df[df['qstart']==1.0] 
	BR1080f = df[df['qstart']==51.0]
	p1_amps = []
	if len(BMEI1427) >= 1 and len(BR1080f) >= 1:
		for i, BME_row in BMEI1427.iterrows():
			BME_start = BME_row['sstart']
			BME_end = BME_row['send']
			for j, BR_row in BR1080f.iterrows():
				BR_start = BR_row['sstart']
				BR_end = BR_row['send']
				amp_start, amp_end = min_max(BME_start, BME_end, BR_start, BR_end)
				amplicon = amp_end - amp_start
				p1_amps.append(amplicon)
		return(p1_amps)
	else:
		return(np.inf)

def primer_two(df, sseqid):
	BR1080r = df[df['qstart']==1.0]
	BMEI1688 = df[df['qstart']==51.0]
	#print(sseqid)
	#print(BR1080r)
	#print('\n')
	#print(BMEI1688)
	#print('\n')
	p2_amps = []
	if len(BR1080r) >= 1 and len(BMEI1688) >=1:
		for i, BR_row in BR1080r.iterrows():
			BR_start = BR_row['sstart']
			BR_end = BR_row['send']
			for j, BME_row in BMEI1688.iterrows():
				BME_start = BME_row['sstart']
				BME_end = BME_row['send']
				amp_start, amp_end = min_max(BME_start, BME_end, BR_start, BR_end)
				amplicon = amp_end - amp_start
				p2_amps.append(amplicon)
		return(p2_amps)
	else:
		return(np.inf)
'''

def primer_one(df, sseqid):
	BMEI1426 = df[df['qstart']==1.0] 
	start_26 = BMEI1426['sstart']
	end_26 = BMEI1427['send']
	BMEI1427 = df[df['qstart']==51.0]
	start_27 = BMEI1426['sstart']
	end_27 = BMEI1427['send']
	amp_start, amp_end = min_max(start_26, end_26, start_27, end_27)
	amplicon = amp_end - amp_start
	return(amplicon)

def primer_two(df, sseqid):
	#multicases!! 
	x=0
	
def primer_three(df, sseqid):
	BMEI1688 = df[df['qstart']==1.0] 
	start_88 = BMEI1688['sstart']
	end_88 = BMEI1687['send']
	BMEI1687 = df[df['qstart']==51.0]
	start_87 = BMEI1688['sstart']
	end_87 = BMEI1687['send']
	amp_start, amp_end = min_max(start_88, end_88, start_87, end_87)
	amplicon = amp_end - amp_start	
	return(amplicon)

def primer_four(df, sseqid):
	BMEI0205f = df[df['qstart']==1.0] 
	start_f = BMEI0205f['sstart']
	end_f = BMEI0205r['send']
	BMEI0205r = df[df['qstart']==51.0]
	start_r = BMEI0205f['sstart']
	end_r = BMEI0205r['send']
	amp_start, amp_end = min_max(start_f, end_f, start_r, end_r)
	amplicon = amp_end - amp_start
	return(amplicon)

def idk(sseqid):
	unique_df = blast_df[blast_df['sseqid']==sseqid] # df containing only hits from a single contig
	primer_one_df = unique_df[unique_df['qseqid']=='BMEI1426-BMEI1427'] # subset of unique_df with only primer 1 data
	primer_two_df = unique_df[unique_df['qseqid']=='BR1080f-BR1080r'] # subset of unique_df with only primer 2 data
	primer_three_df = unique_df[unique_df['qseqid']=='BMEI1688-BMEI1687'] # subset of unique_df with only primer 3 data
	primer_four_df = unique_df[unique_df['qseqid']=='BMEI0205f-BMEI0205r'] # subset of unique_df with only primer 4 data
	#return(sseqid, len(primer_one_df), len(primer_two_df), len(primer_three_df), len(primer_four_df))

	if len(primer_one_df) >1:
		p1_amps = primer_one(primer_one_df, sseqid)
	else:
		p1_amps = np.inf

	#if len(primer_two_df)>1:
	#	p2_amps = primer_two(primer_two_df, sseqid)
	#else:
	#	p2_amps = np.inf
	
	if len(primer_three_df) >1:
		p3_amps = primer_three(primer_three_df, sseqid)
	else:
		p3_amps = np.inf

	if len(primer_four_df)>1:
		p4_amps = primer_four(primer_four_df, sseqid)
	else:
		p4_amps = np.inf	

	return(sseqid, p1_amps, p3_amps, p4_amps)
	

df = pd.DataFrame(index = ['BMEI1427-BR1080f', 'BR1080r-BMEI1688'])
for seqid in sseqids:
	sseqid, L1, L3, L4 = idk(seqid)
	print(sseqid)
	print('Primer 1 Amp: ', L1)
	print('Primer 3 Amp: ', L3)
	print('Primer 4 Amp: ', L4)
	print('\n')


#df = df.dropna(axis='columns', how='all')
#print(df)
#df.to_csv('test.csv')




'''
sseqids = list(blast_df['sseqid'].unique())
pairs = {'BMEI1427':'BR1080f','BR1080r':'BMEI1688','BR1080f':'BMEI1427','BMEI1688': 'BR1080r'}

def amplicon_df(df, qseqid, pair):
	for i in range(0,len(pair_df)): # for each entry in pair_df 
		row_pair = pair_df.iloc[i]
		start_pair = row_pair['sstart'] #the start location of the pair in the sample
		amplicon = end_sseqid - start_pair # location between the primer and its pair
		return(qseqid, pair, amplicon)

def amplicon_series(df, qseqid, pair):
	start_pair = pair_df['sstart']
	amplicon = end_sseqid - start_pair # location between the primer and its pair
	return(qseqid, pair, amplicon)

def idk(sseqid):
	unique_df = blast_df[blast_df['sseqid']==sseqid]
	if len(unique_df) > 1: # as long as there is more than one hit in the given sample
		qseqids = list(set(unique_df.index.values)) # a set of the primer names that appear in the sample
		
		for i in range(0,len(unique_df)): # for each row in unique_df
			row_sseqid = unique_df.iloc[i]
			end_sseqid = row_sseqid['send'] # end location provided by blast
			pair = pairs[row_sseqid.name] # returns the primers pair
			#print(row_sseqid.name)
			#print(pair)
			#print('\n')
			

			if pair in qseqids: #as long as that pair is found in the df (aka found in the set generated above)
				pair_df = unique_df.T[pair] # df with only the pair primers
				pair_df = pair_df.T
				if isinstance(pair_df, pd.DataFrame):
					amplicon_df(pair_df, )
					print(pair_df)

	#else:
	#	return('No Pairs Found in Sequence ', sseqid)

idk('000236255_NC_016797.1')'''
