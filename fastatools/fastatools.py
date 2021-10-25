import pandas as pd
import json
from pathlib import Path
from Bio import Entrez
from datasetmaker.datasetmaker import DatasetMaker
from abc import abstractmethod


class FastaDatasetMaker(DatasetMaker):

    def __init__(self):
        self.column_names = ['pdb_code', 'chain',
                             'chain_type', 'description', 'sequence']

    @abstractmethod
    def assign_chain_type(self, df):
        pass

    def transform_dataframe(self, df):
        return (df
                .pipe(self.process_header)
                .pipe(self.collapse_rows, 'chain')
                .pipe(self.assign_chain_type)
                )

    def organize_dataframe(self, df):
        return (df[self.column_names]
                .sort_values(['pdb_code', 'chain'], ignore_index=True))

    def process_header(self, df):
        header = df['header']
        header = header.str.split('|')
        df['pdb_code'] = header.str[1]
        description = header.str[2]
        df['chain'] = description.str[0]
        df['description'] = description.str.split(', ').str[-1]
        return df.drop('header', axis=1)

    def sort_then_join(self, lst):
        return ','.join(sorted(lst))

    # collapse rows that differ only by specified column
    def collapse_rows(self, df, column_to_collapse):
        columns = [n for n in df.columns if n != column_to_collapse]
        return (df
                .groupby(columns)[column_to_collapse]
                .apply(self.sort_then_join)
                .reset_index())


def read_text_file(filepath):
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


# returns a dataframe of sequences and headers from a FASTA file
def make_dataframe_from_fasta(fasta):
    return pd.DataFrame(process_fasta(fasta), columns=['header', 'sequence'])


def clean_dataframe_header(df, replacements):
    df.columns = df.columns.str.strip().str.lower()
    for key, value in replacements.items():
        df.columns = df.columns.str.replace(key, value, regex=False)
    return df


def make_query(pdb_code):
    return f'{pdb_code}[All Fields] AND pdb[filter]'


def make_fasta_id_list(pdb_codes):
    id_list = []
    n_codes = len(pdb_codes)
    for i, pdb_code in enumerate(pdb_codes):
        query = make_query(pdb_code)
        search_handle = Entrez.esearch(
            db="protein", term=query, idtype="acc", usehistory="y", retmax=50)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        for _id in search_results['IdList']:
            id_list.append(_id)
        print(f'PDB code: {pdb_code}: {i+1} / {n_codes}')
    print(f'Search for {n_codes} proteins is complete')
    return id_list


def get_fasta_id_list(pdb_codes):
    print(f'Searching for {len(pdb_codes)} PDB entries')
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
    id_list = get_fasta_id_list(pdb_codes)
    print(f'Number of IDs (chains) found: {len(id_list)}')
    search_handle = Entrez.epost(db='protein', id=','.join(map(str, id_list)))
    search_results = Entrez.read(search_handle, validate=True)
    search_handle.close()
    return search_results, id_list


def fetch_fasta(search_results, id_count):
    fetch_handle = Entrez.efetch(
        db="protein",
        rettype="fasta",
        retmode="text",
        retstart=0,
        retmax=id_count,
        webenv=search_results["WebEnv"],
        query_key=search_results["QueryKey"],
        idtype="acc"
    )
    data = fetch_handle.read()
    fetch_handle.close()
    return data


def get_fasta_from_ncbi_query(pdb_codes, email, api_key):
    Entrez.email = email
    Entrez.api_key = api_key
    search_results, id_list = get_ncbi_search_results(pdb_codes)
    id_count = len(id_list)
    print(f'Downloading {id_count} records')
    return fetch_fasta(search_results, id_count)


def write_fasta_from_ncbi_query(filepath, pdb_codes, email, api_key):
    data = get_fasta_from_ncbi_query(pdb_codes, email, api_key)
    with open(filepath, "w") as out_handle:
        out_handle.write(data)
