import pandas as pd
import mysql.connector as mysql
import argparse
import sys

def read_credentials(file:str) -> dict:

    credentials = {}
    with open(file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if len(line.strip()) > 0 and line[0] != '#' and '=' in line:
            entry = line.strip().split('=')
            credentials[entry[0]] = entry[1]

    return credentials

def get_data(db_config:dict, query:str, columns=None, params=None) -> pd.DataFrame:

    # MySQL connection
    try:
        conn = mysql.connect(**db_config)
        cursor = conn.cursor()

        # Execute the query
        if params:
            cursor.execute(query, params)
        else:
            cursor.execute(query)

        myresults = cursor.fetchall()

        if columns:
            myresults=pd.DataFrame(myresults, columns=columns)

        return myresults

    except mysql.Error as err:
        log_msg = f"{err}"
        print(log_msg)
        sys.exit(2)

def get_tables(config, table_name, file_name, type:str='dna') -> pd.DataFrame:

    query = f'describe {type.lower()}_{table_name};'
    table = get_data(config, query)
    columns = [x[0] for x in table]
    
    query = f'select * from {type.lower()}_{table_name};'
    table = get_data(config, query, columns)

    table.to_csv(file_name, sep='\t', index=False)

def retrieve_data(config, prefix):

    for type in ['aa', 'dna']:
        for table in ['alignments', 'alignments_taxonomy']:
            file_name = f'{prefix}{type}_{table}.tsv'
            get_tables(config, table, file_name, type=type)

def main():

    parser = argparse.ArgumentParser(description="Get EvoNAPS tables.")
    parser.add_argument("--config", type=str, required=True, help="Path to the config file holding EvoNAPS database credentials.")
    parser.add_argument("--prefix", type=str, required=False, default='./', help="Option to declare prefix for output file. Default is current directory.")
    args = parser.parse_args()

    if args.prefix[-1] != '/':
        args.prefix += '/'

    credentials = read_credentials(args.config)
    retrieve_data(credentials, args.prefix)       
    
    return 0

if __name__ == "__main__":
    main()