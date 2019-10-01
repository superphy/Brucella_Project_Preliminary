import pandas as pd 
import os 
import ete3
i=0
j=0
ts = ete3.TreeStyle()

metadata = pd.read_csv('../Metadata.csv', dtype="str")
meta_dict = {}
strains = ['B. abortus','B. canis','B. ceti', 'B. inopinata', 'B. melitensis', 'B. microti', 'B. neotomae', 'B. ovis', 'B. pinnipedialis', 'B. sp', 'B. suis', 'no strain found']
col_dict = {'B. abortus': 'OliveDrab','B. canis': 'SaddleBrown','B. ceti':'MediumTurquoise', 'B. inopinata':'Purple', 'B. melitensis':'FireBrick', 'B. microti':'DarkOrange', 'B. neotomae':'Navy', 'B. ovis':'Plum', 'B. pinnipedialis':'SlateGray', 'B. sp':'Gold', 'B. suis':'HotPink', 'no strain found': 'Lime' }
#Creating a dictionary where the key is the tree id and the value is the strain 
while i< len(metadata):
	met = metadata.iloc[i]
	key = met['Tree ID']
	value = met['Strain']
	meta_dict[key]=value
	i=i+1

#cwd = os.getcwd()
#tree_location = cwd +"/kSNP3_Output/tree_tipAlleleCounts.parsimony.tre"
tree = ete3.Tree("/home/ashlynn/Desktop/Fall_2019/Brucella/kSNP3_Output/tree_tipAlleleCounts.parsimony.tre")



#adding the strain as a feature and color coding
for node in tree.traverse():
	node.add_feature('strain', '')
	if node.is_leaf():
		node.name = node.name[0:9]
		strain = meta_dict[node.name]
		node.strain = strain
		if node.name == "000054005":
			node.strain = "B. abortus"
		if strain == "no strain found":
			node.dist = 0.3

tbd = ["000250835","003966035","000292125","000292205","000250775", "000292045"]
for name in tbd:
	node = tree.search_nodes(name=name)
	node[0].detach()

def remove_Bsp():
	# any node whose strain is B. sp is changed to the strain of its sister node
	for node in tree.traverse(strategy = 'postorder'):
			if node.strain == "B. sp":
				sis = node.get_sisters()
				kiddos = node.get_descendants()
				if sis[0].strain != '':
					node.strain = sis[0].strain
					for kid in kiddos:
						kid.strain = sis[0].strain

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

remove_Bsp()
strain_id()
remove_Bsp()
strain_id()
remove_Bsp()
strain_id()

#adding the approporiate color to each node			
for node in tree.traverse():
	if node.strain !='':
		ns = ete3.NodeStyle() 
		ns['bgcolor'] = col_dict[node.strain]
		node.set_style(ns)

#building the legend
while j < len(col_dict):
	ts.legend.add_face(ete3.TextFace(strains[j]+"		", fgcolor = col_dict[strains[j]], fsize = 206), column = j)
	j=j+1

ts.mode = "c"
ts.legend_position = 2
ts.show_branch_length = True
ts.arc_start = 180
ts.arc_span = 180
#tree.render('tree_test.pdf', tree_style=ts)
tree.show(tree_style = ts)
#tree.write(outfile = "new_tree.nw")
