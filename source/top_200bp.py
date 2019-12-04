import pandas as pd 
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
pd.options.mode.use_inf_as_na = True

ranks = pd.read_csv('../Ranks.csv', index_col = 0)
k_loc = pd.read_csv('../Kmer_Locations.csv', index_col = 0)
strain_occ = pd.read_csv('../Strain_Occurances.csv', index_col = 0)
strains = strain_occ.index.values
chroms = {'C1':'C2', 'C2':'C1'}

def rolling(rank_df):
	roll_mean = rank_df.rolling(window=200).mean()
	return(list(roll_mean))

def strain_df_f(strain):
	strain_df = pd.DataFrame(index = ranks.index.values)
	strain_df['Rank'] = ranks[strain+' RANK']
	strain_df['C1 Loc'] = k_loc['C1 '+strain]
	strain_df['C2 Loc'] = k_loc['C2 '+strain]
	strain_df = strain_df.dropna(axis = 0, thresh=2)
	return(strain_df)

def loc_dict_f(strain_df):
	loc_dict = {}
	locs = list(strain_df['C1 Loc'].dropna())+list(strain_df['C2 Loc'].dropna())
	for i in range(0, int(max(locs))):
		loc_dict[i] = (0, np.inf)
	return(loc_dict)

def top_row(strain, chrom):
	strain_df= strain_df_f(strain)
	loc_dict = loc_dict_f(strain_df)
	 
	strain_df = strain_df.drop(columns=chroms[chrom]+' Loc').dropna(axis =0, how='any')
	for index, row in strain_df.iterrows():
		loc_dict[row[chrom+' Loc']]=(row['Rank'], row.name)
	df = pd.DataFrame()
	
	kmers = [value[1] for value in loc_dict.values()]
	ranks = [value[0] for value in loc_dict.values()]
	
	df['Location'] = loc_dict.keys()
	df['Kmer'] = kmers
	df['Ranks'] = ranks

	roll_mean= rolling(df['Ranks'])
	df['Rolling Avg'] = roll_mean
	df = df.sort_values(by=['Rolling Avg'], ascending=False).head(1)
	return(df)

def intermediate(strain):
	c1_top_row = top_row(strain, 'C1')
	c2_top_row = top_row(strain, 'C2')

	max_ra = max(float(c1_top_row['Rolling Avg']), float(c2_top_row['Rolling Avg']))
	if max_ra == float(c1_top_row['Rolling Avg']):
		out_col = [c1_top_row['Kmer'].to_string(index = False), 'C1',int(c1_top_row['Location'])-200, int(c1_top_row['Location']), float(c1_top_row['Ranks']), float(c1_top_row['Rolling Avg'])]
		return(strain, out_col)
	if max_ra == float(c2_top_row['Rolling Avg']):
		out_col = [c2_top_row['Kmer'].to_string(index = False), 'C2', int(c2_top_row['Location'])-200, int(c2_top_row['Location']), float(c2_top_row['Ranks']), float(c2_top_row['Rolling Avg'])]
		return(strain, out_col)

final_df = pd.DataFrame(index = ['Kmer', 'Chromosome','Start Location', 'End Location', 'Rank', 'Rolling Avg'])
with ProcessPoolExecutor(cpu_count()-10) as ppe:
	for strain, column in ppe.map(intermediate, strains):
		final_df[strain] = column

final_df.to_csv('Top_200_bp.csv')