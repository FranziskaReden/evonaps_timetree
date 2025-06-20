from Bio import Phylo
import json
import pandas as pd
import logging

class TimeTree:

    def __init__(self, folder):
        
        self.tree, self.leaves = None, None
        self.lineages_dict = {}
        self.reversed_dict = {}
        self.folder = folder

    def read_tax_ids(self, file_name:str) :
        tax_ids = pd.read_csv(file_name, sep="\t")
        tax_ids_dict = {}
        for _, row in tax_ids.iterrows():
            tax_ids_dict[row['name']] = row['tax_id']

        return tax_ids_dict

    def rename_tree(self, taxa_file:str) -> None:

        tax_ids = self.read_tax_ids(taxa_file)

        failed = []
        for clade in self.leaves:
            if clade.name:
                if not clade.name.startswith("*"):
                    if str(clade.name) in tax_ids.keys():
                        clade.name = str(tax_ids[str(clade.name)])
                    else:
                        print(f'Warning: taxid for {clade.name} was not found.')

            else: 
                failed.append(clade)
                print(clade)

    def write_tree(self, target_file:str):
        Phylo.write(self.tree, f'{self.folder}{target_file}', format="newick")

    def read_tree(self, tree_file:str):
        self.tree = Phylo.read(tree_file, 'newick', rooted=True)

    def read_lineages(self, lineages_file:str, reversed_ineages_file:str):
        self.lineages_dict = read_json(lineages_file)
        self.reversed_dict = read_json(reversed_ineages_file)

    def set_inner_node(self, tax_id:str, mrca):
        node_name = mrca.name
        mrca.name = tax_id
        self.lineages_dict[tax_id]['included'] = 1
        self.get_age(mrca.name)
        logging.info(f"Setting inner node {node_name} to ID {tax_id} with age {self.lineages_dict[tax_id]['age']}")

    def replace_neighbour(self, tax_id:str, child:str, ancestor:str) -> None:
        self.lineages_dict[ancestor]['neighbours'] = [child if x == tax_id else x for x in self.lineages_dict[ancestor]['neighbours']]
        self.lineages_dict[tax_id]['merged'].append(ancestor)
        self.lineages_dict[tax_id]['included'] = -1
        logging.info(f"Merging {tax_id} with {ancestor}")

    def remove_from_neighbours(self, id:str, ancestor:str) -> None:
        self.lineages_dict[ancestor]['neighbours'] = [item for item in self.lineages_dict[ancestor]['neighbours'] if item != id]

    def combine_neighbours(self, id:str, ancestor:str):
        self.lineages_dict[ancestor]['neighbours'] += self.lineages_dict[id]['neighbours']

    def clean_up_lineages(self, taxa:list, ancestor:str) -> None:

        for lin in taxa:
            self.lineages_dict[lin]['included'] = -1
            self.lineages_dict[lin]['merged'].append(ancestor)
            logging.info(f"Merging {lin} with {ancestor}")

        # Remove last entry in the list (direct child) from ancestor neighbour
        self.remove_from_neighbours(taxa[-1], ancestor)
        # Instead add neighbours (children) of the tax_id in question (all resolved) to the neighbours of the ancestor
        self.combine_neighbours(taxa[0], ancestor)

    def get_lineage_back(self, id):
        lineage = [id]
        clade = lineage[0]

        while clade != '1':
            clade = self.reversed_dict[clade]
            lineage.append(clade)

        return lineage

    def check_ancestry(self, tax1:str, tax2:str) -> tuple[str, list, list]:

        lineage1 = self.get_lineage_back(tax1)
        lineage2 = self.get_lineage_back(tax2)

        for clade in lineage1:
            if clade in lineage2:
                mrca = clade
                break

        return mrca, lineage1[:lineage1.index(mrca)], lineage2[:lineage2.index(mrca)]
    
    def get_leaf(self, id):

        candidate = id
        if self.lineages_dict[candidate]['leaf'] == 1:
            return None
        
        while self.lineages_dict[candidate]['leaf'] != 1:
            candidate = self.lineages_dict[candidate]['neighbours'][0]
        
        return candidate

    def get_distance(self, inner_node, leaf):       
        return self.tree.distance(inner_node, leaf)
    
    def get_age(self, inner_node):

        leaf = self.get_leaf(inner_node)
        if leaf is None:
            return
        
        self.lineages_dict[inner_node]['age'] = self.get_distance(inner_node, leaf)

    def check_resolved_taxa(self) -> list:
        
        to_be_resolved = []

        for tax_id in self.lineages_dict.keys():
            if self.lineages_dict[tax_id]['included'] == 0:
                in_lineage = self.lineages_dict[tax_id]['neighbours']
                tmp = True
                for taxon in in_lineage:
                    if self.lineages_dict[str(taxon)]['included'] == 0:
                        tmp = False
                        #print(f'for tax_id {tax_id}, tax_id {taxon} is not resolved.')
                        break
                if tmp is True:
                    to_be_resolved.append(tax_id)

        return to_be_resolved

def read_json(file:str) -> dict:
    """
    Reads a JSON file and returns a dictionary.
    """
    with open(file, 'r') as f:
        return json.load(f)

def write_json(dict:dict, file:str) -> None:
    """
    Writes a dictionary to a JSON file.
    """
    with open(file, 'w') as f:
        json.dump(dict, f, indent=4)

def read_lineage(line:str) -> list:
    """
    Parses a lineage string and returns a dictionary of tax IDs and names.
    The input string is expected to be in the format "tax_id:name;tax_id:name; ...".
    """
    tax_ids = []
    ranks = line.split(";")
    if ranks[-1] == "":
        ranks = ranks[:-1] # Remove the last empty element
    for rank in ranks:
        tax_id = rank.split(":")[-1]
        tax_ids.append([rank.replace(f':{tax_id}', ''.strip()), tax_id.strip()]) 
    
    tax_ids.reverse()
    return tax_ids

def tabulate_names(tree:Phylo.BaseTree.Tree, len_external:int) -> tuple[dict, list]:
    """
    Tabulates names from the tree, assigning unique names to internal nodes and cleaning up leaf names.
    Returns a dictionary mapping names to clades and a list of names.
    """

    internal_index = len_external + 1
    for clade in tree.find_clades(): 
        if not clade.name: 
            clade.name = "*"+str(internal_index)
            internal_index += 1

def read_timetree(tree_file) -> tuple[Phylo.BaseTree.Tree, list]:
    """	
    Reads a tree file in Newick format and extracts leaf names.
    Returns a dictionary of names (internal and external), a list of leave names, and the tree object.
    """

    tree = Phylo.read(tree_file, "newick", rooted=True)
    external = tree.get_terminals() # Get all leafes (terminals).
    tabulate_names(tree, len(external))
    leaf_names = tree.get_terminals()
    
    return tree, leaf_names