import pandas as pd
import argparse
from tqdm import tqdm
from Bio import Phylo
import logging

from utils import *

def initialize(args, folder):

    data = TimeTree(folder)
    data.tree, data.leaves = read_timetree(args.tree)
    data.rename_tree(args.tax_ids)
    data.write_tree(f'TimeTree5_renamed.nwk')

    data.read_lineages(args.lineages, args.lineages.replace('.json', '_reversed.json'))

    logging.basicConfig(
        filename=f'{folder}retrieve_age.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    return data

def resolve_taxa(data:TimeTree) -> None:

    while True:
        if not get_mrca(data):
            break
        
        # Add breakpoint, write temporary results to files
        write_json(data.lineages_dict, f'{data.folder}TimeTree5_lineages_resolved.json')
        Phylo.write(data.tree, f'{data.folder}TimeTree5_renamed_resolved.nwk', format="newick")

    # Finally, set root if not in tree
    if data.lineages_dict['1']['resolved'] == -1:
        neighbours = data.lineages_dict['1']['neighbours']
        if len(neighbours) == 1 and data.lineages_dict[neighbours[0]]['resolved'] == 1:
            tmp_clade = data.tree.find_any(neighbours[0])
            data.set_inner_node('1', tmp_clade)
            data.lineages_dict[neighbours[0]]['resolved'] = -1

def get_mrca(data:TimeTree) -> bool:

    # Check all clades that can be resolved
    to_be_resolved = data.check_resolved_taxa()
    if to_be_resolved == []:
        return False
    
    for tax_id in tqdm(to_be_resolved, total=len(to_be_resolved), desc="Resolving MRCA"):
        if data.lineages_dict[tax_id]['included'] == 0:
            check_neighbours(tax_id, data)

    return True
    
def check_neighbours(tax_id:str, data:TimeTree):

    # Get the children names (of leafes or inner nodes) in the subtree
    names = data.lineages_dict[tax_id]['neighbours']

    # If there is only one name (no bifurication in tree)
    if len(names) == 1:
        ancestor = data.reversed_dict[tax_id]
        # Replace the tax_id with the child name in the neighbours list of the ancestor
        data.replace_neighbour(tax_id, names[0], ancestor)

    else:
        # Find inner node that is the most recent common ancestor (MRCA) of the names
        mrca = data.tree.common_ancestor(names)

        # If found node is not resolved yet, assign tax_id to it
        if mrca.name.startswith("*"):
            data.set_inner_node(tax_id, mrca)
        
        # If the inner node is already resolved...
        elif mrca.name.startswith("_"):
            logging.warning(f"Unable to resolve {mrca.name} with {tax_id}")
        
        else:
            resolve_inner_node(tax_id, mrca, data)

def resolve_inner_node(tax_id, mrca, data:TimeTree):

    # Check each lineage and find common ancestor
    ancestor, lin1, lin2 = data.check_ancestry(mrca.name, tax_id)
    
    if ancestor != '1':
        # Write out conflicting tax ids and their ancestor
        logging.info(f"Iner node {mrca.name} conflicts with ID {tax_id} - will be set to {ancestor}")

        # Check if found node is actually common ancestor.
        if ancestor == mrca.name:
            # Then just clean up other lineage
            data.clean_up_lineages(lin2, ancestor)

        # Otherwise, clean up both lineages and set found node to ancestor
        else:
            if ancestor == tax_id:
                data.clean_up_lineages(lin1, ancestor)
            
            else:
                if len(lin1) == 0 or len(lin2) == 0:
                    logging.warning(f'Warning: lineages are of size zero for {mrca.name} and {tax_id} with ancestor {ancestor}')
                    return
                
                data.clean_up_lineages(lin1, ancestor)
                data.clean_up_lineages(lin2, ancestor)

            if data.lineages_dict[ancestor]['included'] == 1:
                logging.info(f"ID {ancestor} already in tree. Will be set to current node {mrca.name} instead")
                # Rename clade with ancestor names
                while (tmp_clade := data.tree.find_any(ancestor)):
                    tmp_clade.name = f'**{ancestor}'

            # Finally, set found node to ancestor.
            data.set_inner_node(ancestor, mrca)
    
    elif tax_id == '1':
        data.clean_up_lineages(lin1, ancestor)
        data.set_inner_node(tax_id, mrca)

    elif mrca.name == '1':
        data.clean_up_lineages(lin2, ancestor)

    # If no common ancestor is found, flag a warning, write out conflicting tax_ids
    else:
        logging.warning(f"Iner node {mrca.name} in conflict with ID {tax_id} cannot be resolved (no ancestor found).")

def main():
    parser = argparse.ArgumentParser(description="Get the most recent common ancestor (MRCA) for each taxon in the TimeTree.")
    parser.add_argument("--lineages", type=str, required=True, help="Path to the lineages JSON file.")
    parser.add_argument("--tree", type=str, required=True, help="Path to the tree file in Newick format.")
    parser.add_argument("--tax_ids", type=str, required=True, help="Path to the tax IDs file in TSV format.")
    parser.add_argument("--prefix", type=str, help="Option to declare output folder.")
    args = parser.parse_args()
    
    folder = './'
    if args.prefix:
        folder = args.prefix
        if not folder.endswith('/'):
            folder += '/'

    data = initialize(args, folder)

    resolve_taxa(data)

    return 0

if __name__ == "__main__":
    main()