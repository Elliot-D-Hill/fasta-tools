import pandas as pd
import json
from pathlib import Path
from Bio import Entrez
from abc import ABC, abstractmethod


class FastaDatasetMaker(ABC):
    def __init__(self, raw_path, cleaned_path, fasta_path, processed_path, email, ncbi_api_key):
        self.raw_path = raw_path
        self.cleaned_path = cleaned_path
        self.fasta_path = fasta_path
        self.processed_path = processed_path
        self.email = email
        self.ncbi_api_key = ncbi_api_key
        self.column_names = ['pdb_code', 'chain',
                             'chain_type', 'description', 'sequence']
        self.make_cleaned_dataframe()
        self.df = read_dataframe(self.cleaned_path)
        self.pdb_codes = self.get_pdb_codes()

    @abstractmethod
    def filter_dataframe(self, df):
        pass

    @abstractmethod
    def format_dataframe(self, df):
        pass

    @abstractmethod
    def assign_chain_type(self, df):
        pass

    def make_cleaned_dataframe(self):
        (read_dataframe(self.raw_path)
            .pipe(self.format_dataframe)
            .pipe(self.filter_dataframe)
            .reset_index(drop=True)
            .to_csv(self.cleaned_path, index=False)
         )

    def get_pdb_codes(self):
        return self.df['pdb_code'].str[0:4].unique()

    def make_dataset(self):
        write_fasta_from_ncbi_query(
            self.pdb_codes, self.fasta_path, self.email, self.ncbi_api_key)
        (make_fasta_dataframe(self.fasta_path)
            .pipe(process_header)
            .pipe(self.assign_chain_type)
            .pipe(collapse_rows, 'chain')
            .pipe(organize_fasta_dataframe, self.column_names)
            .to_csv(self.processed_path, index=False)
         )


def read_dataframe(filepath):
    extension = str(filepath).split('.')[-1]
    if extension == 'tsv':
        return pd.read_csv(filepath, sep='\t')
    elif extension == 'csv':
        return pd.read_csv(filepath, sep=',')


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
def make_fasta_dataframe(filepath):
    fasta = read_text_file(filepath)
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
        print(f'PDB code: {pdb_code}: {i+1} / {n_codes}')
        search_handle = Entrez.esearch(
            db="protein", term=query, idtype="acc", usehistory="y", retmax=50)
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


def write_fasta_from_ncbi_fetch(search_results, id_list, filepath):
    id_count = len(id_list)
    print(f'Downloading {id_count} records')
    with open(filepath, "w") as out_handle:
        data = fetch_fasta(search_results, id_count)
        out_handle.write(data)


def write_fasta_from_ncbi_query(pdb_codes, filepath, email, api_key):
    Entrez.email = email
    Entrez.api_key = api_key
    search_results, id_list = get_ncbi_search_results(pdb_codes)
    write_fasta_from_ncbi_fetch(search_results, id_list, filepath)


def process_header(df):
    header = df['header']
    header = header.str.split('|')
    df['pdb_code'] = header.str[1]
    description = header.str[2]
    df['chain'] = description.str[0]
    df['description'] = description.str.split(', ').str[-1]
    return df[['pdb_code', 'chain', 'description', 'sequence']]


def organize_fasta_dataframe(df, column_names):
    df = df[column_names]
    return df.sort_values(['pdb_code', 'chain'], ignore_index=True)


def sort_then_join(lst):
    return ','.join(sorted(lst))


# collapse rows that differ only by one column
def collapse_rows(df, column_to_collapse):
    columns = [n for n in df.columns if n != column_to_collapse]
    df = df.groupby(columns)[column_to_collapse].apply(sort_then_join)
    return df.reset_index()
