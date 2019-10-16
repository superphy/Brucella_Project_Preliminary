import pandas as pd 
import os 
import ete3
i=0
j=0
k = 1
ts = ete3.TreeStyle()

metadata = pd.read_csv('Metadata.csv', dtype="str")
meta_dict = {}
strains = ['Brucella suis','Brucella canis','Brucella sp','Brucella neotomae','Brucella pinnipedialis','Brucella ceti', 'Brucella ovis','Brucella abortus','Brucella melitensis']
col_dict = {'Brucella sp':'#ebc04b', 'Brucella suis':'#f8a4d0','Brucella abortus': '#72a24a','Brucella canis': '#b8642f','Brucella ceti':'#59cf9e', 'Brucella melitensis':'#e8382e', 'Brucella neotomae':'#6c86e2', 'Brucella ovis':'#ba6eb4', 'Brucella pinnipedialis':'#f4792a' }

#Creating a dictionary where the key is the tree id and the value is the strain 
while i< len(metadata):
	met = metadata.iloc[i]
	sample = met['Sample']
	key = sample[4:13]
	value = met['Species']
	meta_dict[key]=value
	i=i+1

cwd = os.getcwd()
tree_location = cwd +"/kSNP3_Output/tree.parsimony.tre"
tree = ete3.Tree(tree_location)
#tree = ete3.Tree("/home/ashlynn/Desktop/Fall_2019/Brucella/kSNP3_Output/tree.parsimony.tre")

#adding the strain as a feature and color coding
for node in tree.traverse():
	node.add_feature('strain', '')
	if node.is_leaf():
		strain = meta_dict[node.name]
		node.strain = strain
	else: 
		node.name = k
		k+=1

def strain_id():
	#For every node in the tree that is not a leaf, if all its decendants are the same strain, let that node be of the same strain
	for node in tree.traverse():
		if node.is_leaf() == False:
			decs = node.get_descendants()
			strains = []
			for dec in decs:
				if dec.is_leaf():
					strains.append(dec.strain)
			if all(x==strains[0] for x in strains):
				node.add_feature('strain', strains[0])

strain_id()

#adding the approporiate color to each node			
for node in tree.traverse():
	if node.name == 12:
		node.swap_children()
	if node.name == 20:
		node.swap_children()
	if node.strain !='':
		ns = ete3.NodeStyle()
		ns['bgcolor'] = col_dict[node.strain]
		node.set_style(ns)
	if node.dist > 0.02 :
		node.dist = node.dist*0.01


#building the legend
while j < len(col_dict):
	ts.legend.add_face(ete3.TextFace(strains[j]+"		", fgcolor = col_dict[strains[j]], fsize = 600), column = j)
	j=j+1

ts.mode = "c"
ts.legend_position = 2
ts.show_branch_length = True
ts.arc_start = 180
ts.arc_span = 180
tree.render('tree.pdf', tree_style=ts)
#tree.show(tree_style = ts)

