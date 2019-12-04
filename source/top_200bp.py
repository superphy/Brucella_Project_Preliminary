import pandas as pd 
import numpy as np
pd.options.mode.use_inf_as_na = True

ranks = pd.read_csv('../Ranks.csv', index_col = 0)
k_loc = pd.read_csv('../Kmer_Locations.csv', index_col = 0)
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

def c1(strain, chrom):
	strain_df= strain_df_f(strain)
	loc_dict = loc_dict_f(strain_df)
	 
	strain_df = strain_df.drop(columns=chroms[chrom]+' Loc').dropna(axis =0, how='any')
	for index, row in strain_df.iterrows():
		loc_dict[row[chrom+' Loc']]=(row['Rank'], row.name)
	df = pd.DataFrame(index = loc_dict.keys())
	
	
	print(loc_dict.values())
	#df['Kmer'] = loc_dict.values()[1]
	#df['Location'] = loc_dict.keys()
	#df['Ranks'] = loc_dict.values()[0]

	#roll_mean= rolling(df['Ranks'])
	#df['Rolling Avg'] = roll_mean
	#df = df.sort_values(by=['Rolling Avg'], ascending=False)
	return(strain_df)

def c2(strain):
	strain_df = strain_df_f(strain)
	strain_df = strain_df.drop(columns='C1 Loc').dropna(axis = 0, how = 'any').sort_values(by=['C2 Loc'])
	#
	#roll_mean= rolling(strain_df['Rank'])
	#strain_df['Rolling Avg'] = roll_mean
	#strain_df.to_csv('test.csv')
	return(strain_df)



c1('Brucella abortus', 'C1')
