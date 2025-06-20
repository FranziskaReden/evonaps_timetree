import pandas as pd
import argparse
import sys
import os
from tqdm import tqdm

from utils import read_json

def check_files(files:list):

    for file in files:
        if os.path.isfile(file) is False:
            print(f'File {file} was not found!')
            sys.exit(2)

def check_merged_calde(id:str, ages:dict) -> str:
    
    while ages[id]['included'] != 1:
        if len(ages[id]['merged']) < 1:
            return None
        id = ages[id]['merged'][-1]

    return id

def get_ages(ages:dict, tax_table:pd.DataFrame, file_name:str):
    
    tax_table['AGE'] = None
    tax_table['MERGED'] = None
    ages['1']['age'] = 3772.466

    tax_ids = tax_table['LCA_TAX_ID'].unique()
    
    for tax_id in tqdm(tax_ids, total=len(tax_ids), desc='Getting ages...'):
        if str(tax_id) in ages.keys():
            indeces = tax_table[tax_table['LCA_TAX_ID'] == tax_id].index
            if ages[str(tax_id)]['included'] == 1:
                for idx in indeces:
                    tax_table.at[idx, 'AGE'] = float(ages[str(tax_id)]['age'])

            else:
                merged_id = check_merged_calde(str(tax_id), ages)
                if merged_id is not None:
                    for idx in indeces:
                        tax_table.at[idx, 'AGE'] = float(ages[str(merged_id)]['age'])
                        tax_table.at[idx, 'MERGED'] = int(merged_id)

    tax_table.to_csv(file_name, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Get mrca anges for EvoNAPS alignments.")
    parser.add_argument("--prefix", type=str, required=True, default='./', help="Option to declare prefix for output file. Default is current directory.")
    args = parser.parse_args()

    if args.prefix[-1] != '/':
        args.prefix += '/'
    
    # Declare file names
    ages_file = f'{args.prefix}TimeTree5_lineages_resolved.json'
    aa_ali_file = f'{args.prefix}aa_alignments.tsv'
    aa_tax_file = f'{args.prefix}aa_alignments_taxonomy.tsv'
    dna_ali_file = f'{args.prefix}dna_alignments.tsv'
    dna_tax_file = f'{args.prefix}dna_alignments_taxonomy.tsv'

    # Check if files exist
    check_files([ages_file, aa_ali_file, aa_tax_file, dna_ali_file, dna_tax_file])

    # Read in files
    ages = read_json(ages_file)
    aa_ali = pd.read_csv(aa_ali_file, sep='\t', low_memory=False)
    aa_tax = pd.read_csv(aa_tax_file, sep='\t', low_memory=False)
    dna_ali = pd.read_csv(dna_ali_file, sep='\t', low_memory=False)
    dna_tax = pd.read_csv(dna_tax_file, sep='\t', low_memory=False)

    get_ages(ages, aa_tax, aa_tax_file.replace('.tsv', '_ages.tsv'))
    get_ages(ages, dna_tax, dna_tax_file.replace('.tsv', '_ages.tsv'))


if __name__ == "__main__":
    main()