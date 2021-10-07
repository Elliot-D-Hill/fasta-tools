from Bio import Entrez
import pandas as pd
from pathlib import Path
import json


def read_file(filepath):
    with open(filepath, "r") as f:
        return f.read()


# returns a tuble containing the sequence header and the sequence itself
def separate_header(sequence):
    return sequence[0], ''.join(sequence[1:])


# returns a list of tuples containing a header and sequence for each sequence in a FASTA file
def process_fasta(file_text):
    sequences = file_text.split('>')[1:]
    sequences = [sequence.split('\n') for sequence in sequences]
    return [separate_header(sequence) for sequence in sequences]


# returns lists sequences and headers from FASTA file
def make_fasta_dataframe(filepath):
    fasta = read_file(filepath)
    return pd.DataFrame(process_fasta(fasta), columns=['header', 'sequence'])


def clean_df_columns(df, replacements):
    df.columns = df.columns.str.strip().str.lower()
    for key, value in replacements.items():
        df.columns = df.columns.str.replace(key, value, regex=False)


def make_query(pdb_code):
    return f'{pdb_code}[All Fields] AND pdb[filter]'


def make_fasta_id_list(pdb_codes):
    id_list = []
    n_codes = len(pdb_codes)
    for i, pdb_code in enumerate(pdb_codes):
        query = make_query(pdb_code)
        print(f'PDB code: {pdb_code}: {i+1} / {n_codes}')
        search_handle = Entrez.esearch(
            db="protein", term=query, idtype="acc", usehistory="y")
        search_results = Entrez.read(search_handle)
        search_handle.close()
        for _id in search_results['IdList']:
            id_list.append(_id)
    return id_list


def get_fasta_id_list(pdb_codes):
    id_list_filepath = Path('Data/Raw/fasta_id_list.txt')
    if not id_list_filepath.is_file():
        id_list = make_fasta_id_list(pdb_codes)
        with open(id_list_filepath, 'w') as f:
            f.write(json.dumps(id_list))
    else:
        with open(id_list_filepath, 'r') as f:
            id_list = json.loads(f.read())
    return id_list


def get_ncbi_search_results(pdb_codes):
    print(f'Searching for {len(pdb_codes)} PDB entries')
    id_list = get_fasta_id_list(pdb_codes)
    print(f'Number of IDs (chains): {len(id_list)}')
    search_handle = Entrez.epost(db='protein', id=','.join(map(str, id_list)))
    search_results = Entrez.read(search_handle, validate=True)
    search_handle.close()
    return search_results, id_list


def write_fasta_from_ncbi_query(pdb_codes, filepath):
    Entrez.email = "ehill3@bidmc.harvard.edu"
    Entrez.api_key = 'd8b7d4b6c3abcadbea8bcd405fa992c55409'
    search_results, id_list = get_ncbi_search_results(pdb_codes)
    id_count = len(id_list)
    print(f'Downloading {id_count} records')
    with open(filepath, "w") as out_handle:
        fetch_handle = Entrez.efetch(
            db="protein",
            rettype="fasta",
            retmode="text",
            retstart=0,
            retmax=len(id_list),
            webenv=search_results["WebEnv"],
            query_key=search_results["QueryKey"],
            idtype="acc"
        )
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)


def process_header(df):
    header = df['header']
    header = header.str.split('|')
    df['pdb_code'] = header.str[1]
    df['chain'] = header.str[2].str[0]
    df['description'] = header.str[2].str.split(', ').str[1]
    return df[['pdb_code', 'chain', 'description', 'sequence']]


def save_dataframe(filename, df):
    df.to_csv(filename, index=False)


def sort_then_join(lst):
    return ','.join(sorted(lst))


def reformat_fasta_dataframe(df, column_names):
    df = df[column_names]
    return df.sort_values(['pdb_code', 'chain'], ignore_index=True)


# collapse rows that different only by chains
def collapse_rows(df, row_to_collapse):
    columns = [n for n in df.columns if n != row_to_collapse]
    df = df.groupby(columns)[row_to_collapse].apply(sort_then_join)
    return df.reset_index()
