import skbio.io
import pandas as pd

blast_output_path = '../blast/blast_search_output.tsv'
with open(blast_output_path) as fh:
	blast_df = skbio.io.read(fh, format="blast+6",into=pd.DataFrame,default_columns=True)

blast_df.to_csv('../Blast_Output.csv')
