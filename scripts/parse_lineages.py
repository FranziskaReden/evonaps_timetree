import pandas as pd
import argparse
from tqdm import tqdm

from utils import read_lineage, write_json

def create_lineages_dict(lineages:pd.DataFrame) -> tuple[dict,dict]:

    reverse_dict = {}
    lineages_dict = {}

    for _, row in lineages.iterrows():
        name = row['lineage'].split(';')[-2].split(':')[0]
        lineages_dict[row['tax_id']] = {'included': 1, 'name': name, 'neighbours': [], 'age': 0, 'merged': [], 'leaf': 1}

    for _, row in tqdm(lineages.iterrows(), total=lineages.shape[0], desc="Processing lineages"):
        tmp_list = read_lineage(str(row['lineage']))
        original = str(row['tax_id'])
        for i in range (1, len(tmp_list)):
            tax_id = tmp_list[i][1]
            
            if tax_id not in lineages_dict.keys():
                if i == 1:
                    lineages_dict[tax_id] = {'included': 0, 'name':tmp_list[i][0], \
                                         'neighbours': [original], 'age': None, 'merged': [], 'leaf': 0}
                    reverse_dict[original] = tax_id
                else:
                    lineages_dict[tax_id] = {'included': 0, 'name':tmp_list[i][0], \
                                         'neighbours': [tmp_list[i-1][1]], 'age': None, 'merged': [], 'leaf': 0}
                    reverse_dict[tmp_list[i-1][1]] = tax_id

            elif tmp_list[i-1][1] not in lineages_dict[tax_id]['neighbours']:
                lineages_dict[tax_id]['neighbours'].append(tmp_list[i-1][1])
                lineages_dict[tax_id]['included'] = 0
                reverse_dict[tmp_list[i-1][1]] = tax_id

    reverse_dict['1'] = '1'

    return lineages_dict, reverse_dict

def main():
    parser = argparse.ArgumentParser(description="Read in the lineages and assign tax_ids.")
    parser.add_argument("--lineages", type=str, help="Path to the lineages file to parse.")
    parser.add_argument("--output", type=str, default="data/TimeTree5_lineages_unresolved", \
                        help="Output file for the lineages JSON.")
    args = parser.parse_args()

    print("Reading in the lineages file...")
    lineages = pd.read_csv(args.lineages, sep='\t')

    print("Filtering through lineages...")
    lineages_dict, reverse_dict = create_lineages_dict(lineages)
    write_json(lineages_dict, f'{args.output}.json')
    write_json(reverse_dict, f'{args.output}_reversed.json')
    
    return 0

if __name__ == "__main__":
    main()