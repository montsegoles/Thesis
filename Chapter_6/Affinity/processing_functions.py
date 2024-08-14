from typing import List

import pandas as pd
import numpy as np
import requests
from pandas import DataFrame
from collections import defaultdict

def filter_canonical_and_length_db(max_length: int, df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters dataframe rows based on sequence length and canonical amino acids.
    :param max_length: Maximum length allowed for sequences.
    :param df: DataFrame containing 'protein_1_seq' and 'protein_2_seq' columns.
    :return: Filtered DataFrame.
    """
    def has_noncanonical_db(sequence: str):
        """
        Checks if a protein sequence contains non-canonical amino acids.
        Returns True if non-canonical amino acids are found, False otherwise.
        """
        canonical_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        return any(aa not in canonical_amino_acids for aa in sequence)

    def is_long_sequence_db(sequence, max_length_n):
        """
        Checks if a protein sequence is longer than the specified maximum length.
        Returns True if the sequence is longer, False otherwise.
        """
        return len(sequence) > max_length_n

    # Apply the filtering functions to 'protein_1_seq' and 'protein_2_seq' columns
    df['Has_Noncanonical_Protein_1'] = df['protein_1_seq'].apply(has_noncanonical_db)
    df['Has_Noncanonical_Protein_2'] = df['protein_2_seq'].apply(has_noncanonical_db)
    
    # Create masks for filtering
    mask_protein_1 = ~df['Has_Noncanonical_Protein_1'] & ~df['protein_1_seq'].apply(is_long_sequence_db, args=(max_length,))
    mask_protein_2 = ~df['Has_Noncanonical_Protein_2'] & ~df['protein_2_seq'].apply(is_long_sequence_db, args=(max_length,))
    
    # Apply masks to filter DataFrame
    filtered_df = df[mask_protein_1 & mask_protein_2].copy()
    
    # Drop the temporary columns
    filtered_df.drop(['Has_Noncanonical_Protein_1', 'Has_Noncanonical_Protein_2'], axis=1, inplace=True)
    
    return filtered_df


def has_noncanonical(sequence: str):
    """
    Checks if a protein sequence contains non-canonical amino acids.
    Returns True if non-canonical amino acids are found, False otherwise.
    """
    canonical_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    return any(aa not in canonical_amino_acids for aa in sequence)


def is_long_sequence(sequence, max_length_n):
    """
    Checks if a protein sequence is longer than the specified maximum length.
    Returns True if the sequence is longer, False otherwise.
    """
    return len(sequence) > max_length_n


def filter_canonical_and_length(max_length: int, df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters dataframes according to lengths
    :param max_length: filters sequences by max length
    :param df: dataframe containing the sequences and ids
    :return: dataframe with only canonical amino acids and lengths under the max
    """
    length = max_length
    df['Has_Noncanonical'] = df['protein_seq'].apply(has_noncanonical)
    canonical_df = df.loc[~df['Has_Noncanonical']]
    filtered_df = canonical_df.loc[~canonical_df['protein_seq'].apply(is_long_sequence, max_length_n=length)]
    filtered_df = filtered_df.drop(columns=['Has_Noncanonical'])  # Drop the 'Has_Noncanonical' column
    return filtered_df


def filter_canonical_and_length_gen(max_length: int, df: pd.DataFrame, col_name: str) -> pd.DataFrame:
    """
    Filters dataframes according to lengths
    :param max_length: filters sequences by max length
    :param df: dataframe containing the sequences and ids
    :return: dataframe with only canonical amino acids and lengths under the max
    """
    length = max_length
    df['Has_Noncanonical'] = df[col_name].apply(has_noncanonical)
    canonical_df = df.loc[~df['Has_Noncanonical']]
    filtered_df = canonical_df.loc[~canonical_df[col_name].apply(is_long_sequence, max_length_n=length)]
    filtered_df = filtered_df.drop(columns=['Has_Noncanonical'])  # Drop the 'Has_Noncanonical' column
    return filtered_df


def get_sequence_from_db(db_archive: str, dataframe: pd.DataFrame):
    """
    gets sequences from existing databases
    uses the PDB complex id
    assumes that there are only two proteins
    """
    database = pd.read_csv(db_archive)
    copia = dataframe.copy()
    protein_1_seq = []
    protein_2_seq = []
    for i in range(0, len(dataframe)):
        id_pdb = dataframe['PDB_id'][i]
        result = database[database['PDB_id'] == id_pdb]
        p1seq = str(result['protein_1_seq'].values)
        p2seq = str(result['protein_2_seq'].values)
        # the moved indexes are to cut [" and "]
        protein_1_seq.append(p1seq[2:-2])
        protein_2_seq.append(p2seq[2:-2])
    copia['protein_1_seq'] = protein_1_seq
    copia['protein_2_seq'] = protein_2_seq
    return copia


def get_sequence_from_pdb(pdb_id: str):
    global seq_cache
    if pdb_id in seq_cache:
        return seq_cache[pdb_id]  # Return sequence from cache if available

    # If not in cache, fetch sequence data from the server
    url = 'https://data.rcsb.org/graphql'
    query = '''
    query structure ($id: String!) {
      entry(entry_id:$id){
          rcsb_id
          entry {
              id
          }
          polymer_entities {
              entity_poly {
                  pdbx_seq_one_letter_code_can
              }
              rcsb_polymer_entity_container_identifiers {
                  entity_id
              }
          }
          rcsb_entry_info {
              polymer_entity_count_protein
          }
      }
    }
    '''
    variables = {"id": pdb_id}
    response = requests.post(url, json={"query": query, "variables": variables})

    if response.status_code == 200:
        data = response.json()
        entry = data['data']['entry']
        if entry is None:
            print(f"No data found for PDB ID {pdb_id}")
            return None
        protein_count = entry['rcsb_entry_info']['polymer_entity_count_protein']
        if protein_count == 2:
            sequences = [poly['entity_poly'].get('pdbx_seq_one_letter_code_can', '') for poly in entry.get('polymer_entities', [])]
            if len(sequences) == 2:
                seq_cache[pdb_id] = sequences  # Cache the sequences
                return sequences
            else:
                print(f"Error: Sequences not retrieved for PDB ID {pdb_id}")
                return None
        else:
            print(f"The entry {pdb_id} is not a dimer.")
            return None
    else:
        print(f'Error: {response.status_code}')
        return None


pdb_cache = {}
seq_cache = {}


def clear_cache():
    global pdb_cache
    global seq_cache
    pdb_cache = {}
    seq_cache = {}


def is_dimer(pdb_id):
    global pdb_cache

    if pdb_id in pdb_cache:
        return pdb_cache[pdb_id]
    url = f'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}/fasta?entity=3'
    response = requests.get(url.format(pdb_id), timeout=5000)

    if response.status_code == 200:
        result = True if response.text.strip() == 'null' else False
        pdb_cache[pdb_id] = result
        return result
    else:
        return 'Failure'


def is_dimer_pdb(id):
    global pdb_cache
    url = 'https://data.rcsb.org/graphql'
    query = '''
    query structure ($id: String!) {
      entry(entry_id:$id){
          rcsb_id
          entry {
              id
          }
          polymer_entities {
              entity_poly {
                  pdbx_seq_one_letter_code_can
              }
              rcsb_polymer_entity_container_identifiers {
                  entity_id
              }
          }
          rcsb_entry_info {
              polymer_entity_count_protein
          }
      }
    }
    '''
    variables = {"id": id}
    response = requests.post(url, json={"query": query, "variables": variables})

    if response.status_code == 200:
        data = response.json()
        entry = data['data']['entry']
        if entry is None:
            print(f"No data found for PDB ID {id}")
            return None
        protein_count = entry['rcsb_entry_info']['polymer_entity_count_protein']
        if protein_count == 2:
            result = True
            pdb_cache[id] = result
            return result
        else:
            return False
    else:
        print(f'Error: {response.status_code}')
        return None


def count_macromolecules(pdb_id: str):
    # counts the number of macromolecules in a complex
    url = 'https://data.rcsb.org/graphql'
    query = '''
    query structure ($id: String!) {
      entry(entry_id:$id){
          rcsb_id
          entry {
              id
          }
          polymer_entities {
              entity_poly {
                  pdbx_seq_one_letter_code_can
              }
              rcsb_polymer_entity_container_identifiers {
                  entity_id
              }
          }
          rcsb_entry_info {
              polymer_entity_count_protein
          }
      }
    }
    '''
    variables = {"id": id}
    response = requests.post(url, json={"query": query, "variables": variables})

    if response.status_code == 200:
        data = response.json()
        entry = data['data']['entry']
        if entry is None:
            print(f"No data found for PDB ID {id}")
            return None
        protein_count = entry['rcsb_entry_info']['polymer_entity_count_protein']
        return protein_count
    else:
        print(f'Error: {response.status_code}')
        return None


def drop_duplicates_by_id(dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    eliminates id duplicates, but for just one column
    :rtype: object
    """
    # dataframe -> dataframe
    # uses the PDB ids to eliminate duplicate items
    new_df = dataframe.drop_duplicates(subset=['PDB_id'])
    new_df.reset_index(inplace=True)
    return new_df


def drop_duplicates_by_sequence(dataframe):
    # dataframe -> dataframe
    # isolates the sequence columns and orders each row by alphabetical order
    # returns the dataframe without duplicate interactions
    # resets indexs
    sequences = dataframe.copy()
    # sort the rows
    sequences[['protein_1_seq_sorted', 'protein_2_seq_sorted']] = sequences.apply(
        lambda row: sorted([row['protein_1_seq'], row['protein_2_seq']]), axis=1, result_type='expand')

    # Drop duplicate interaction pairs
    unique = sequences.drop_duplicates(subset=['protein_1_seq_sorted', 'protein_2_seq_sorted'])
    unique.reset_index(drop=True, inplace=True)
    return unique



def eliminate_nulls(dataframe: pd.DataFrame):
    # dataframe -> dataframe
    # eliminates the rows with null values
    copyss = dataframe.copy()
    for i in range(0, len(copyss)):
        if pd.isna(copyss['protein_1_seq'][i]):
            dataframe.drop([i], inplace=True)

        if type(copyss['protein_1_seq'][i]) != str:
            dataframe.drop([i], inplace=True)

        if type(copyss['protein_2_seq'][i]) != str:
            dataframe.drop([i], inplace=True)

        if pd.isna(copyss['Kd [M]'][i]):
            dataframe.drop([i], inplace=True)

    dataframe.reset_index(drop=True, inplace=True)

    return dataframe


def drop_nulls_by(by_what: str, dataframe: pd.DataFrame):
    # str + dataframe -> dataframe
    drop_df = dataframe.dropna(subset=[by_what])
    drop_df.reset_index(drop=True, inplace=True)
    return drop_df


def process_dataframe(df):
    global pdb_cache
    global seq_cache
    pdb_cache = {}  # Initialize cache
    seq_cache = {}  # Initialize sequence cache
    p1_seq = []
    p2_seq = []

    for index, row in df.iterrows():
        pdb_id = row['PDB_id']
        print('Processing:', pdb_id)

        # Check if in cache, else fetch sequences
        sequences = seq_cache.get(pdb_id)
        if sequences is None:
            sequences = get_sequence_from_pdb(pdb_id)

        if sequences:
            p1_seq.append(sequences[0])
            p2_seq.append(sequences[1])
        else:
            print(f"Error: Sequences not retrieved for PDB ID {pdb_id}")
            df.drop(index=index, inplace=True)  # Delete row if sequences not retrieved or not a dimer

    df['protein_1_seq'] = p1_seq
    df['protein_2_seq'] = p2_seq
    return df


def apply_mutations(df: pd.DataFrame, chain_col_1: str, chain_col_2: str, mutation_col: str) -> pd.DataFrame:
    """
    Applies the specified mutations to the sequences of a dataframe
    Adapted for skempi and proximate
    df: pandas dataframe
    The input dataframe should contain at least: PDB_id, Mutations, Chain 1, Chain 2, Prot 1 seq, Prot 2 seq
    mutation format: native aa | chain | position | mutated aa (ex: YA450W)
    chain_col_1 = name of the column with chain 1
    chain_col_2 = name of the column with chain 2
    mutation_col = name of the column with the mutations
    :returns: pd.DataFrame with only the mutated sequences
    """
    # eliminate nulls
    df.dropna(subset=['protein_1_seq', 'protein_2_seq', mutation_col], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Initialize new lists to store mutated sequences
    mutated_seqs_p1 = []
    mutated_seqs_p2 = []
    df_copy = df.copy()

    # Iterate over each row in the dataframe
    for idx, row in df.iterrows():
        # Split mutations by comma if there are multiple mutations
        mutations = row[mutation_col].split(',')

        # Initialize mutated sequences for protein 1 and protein 2
        mutated_seq_p1 = row['protein_1_seq']
        mutated_seq_p2 = row['protein_2_seq']

        # Iterate over mutations
        for mutation in mutations:
            # Parse mutation information
            if mutation[0] == '"':
                mutation = mutation[1:-1]
            length = len(mutation)
            native_aa = mutation[0]
            position = int(mutation[2:(length - 1)])
            mutated_aa = mutation[length-1]
            chain = mutation[1]

            # Update the appropriate protein sequence with the mutation
            if chain in row[chain_col_1]:
                # check length
                if len(mutated_seq_p1) <= position:
                    print(
                        f'Mutation position exceeds protein lenght: protein length = {len(mutated_seq_p1)}, mutation position={position}')
                # check that it is the correct native amino acid in the mutation place
                elif native_aa in mutated_seq_p1[position - 1]:
                    mutated_seq_p1 = mutated_seq_p1[:position - 1] + mutated_aa + mutated_seq_p1[position:]
                # the position is incorrect, delete sequence
                else:
                    print(f'Native aminoacid {native_aa} does not correspond with expected position. '
                          f'There is {mutated_seq_p1[position]} instead.')
            elif chain in row[chain_col_2]:
                # check length
                if len(mutated_seq_p2) <= position:
                    print(
                        f'Mutation position exceeds protein lenght: protein length = {len(mutated_seq_p2)}, mutation position={position}')
                # check that it is the correct native amino acid in the mutation place
                elif native_aa in mutated_seq_p2[position - 1]:
                    mutated_seq_p2 = mutated_seq_p2[:position - 1] + mutated_aa + mutated_seq_p2[position:]
                # the position is incorrect, delete sequence
                else:
                    print(f'Native aminoacid {native_aa} does not correspond with expected position. '
                          f'There is {mutated_seq_p2[position]} instead.')
            else:  # the chain was not found
                print(f"Chain {chain} does not match either interacting chain in row {idx}. Skipping mutation.")

        # Append mutated sequences to the lists
        mutated_seqs_p1.append(mutated_seq_p1)
        mutated_seqs_p2.append(mutated_seq_p2)

    # Add mutated sequences as new columns to the dataframe
    df_copy['mutated_protein_1_seq'] = mutated_seqs_p1
    df_copy['mutated_protein_2_seq'] = mutated_seqs_p2

    # drop failed mutations
    drop_df = df_copy.copy()
    for j in range(0, len(df_copy)):
        # failed mutation
        if df_copy['protein_1_seq'][j] == df_copy['mutated_protein_1_seq'][j] and str(
                df_copy['protein_2_seq'][j]) == str(df_copy['mutated_protein_2_seq'][j]):
            drop_df.drop([j], inplace=True)
        # successful mutation
            # do nothing

    drop_df.reset_index(drop=True, inplace=True)
    drop_df2 = drop_df.dropna(subset=['mutated_protein_1_seq', 'mutated_protein_2_seq'])
    drop_df2.reset_index(drop=True, inplace=True)

    return drop_df2


def process_mutated_db(dataframe: pd.DataFrame) -> list[DataFrame]:
    """
    from a database with wild type sequences and mutated sequences
    creates a dataframe with aditional column "is mutation"
    has unique columns with sequences whether mutated or not | Kd mutated or not
    only keeps Kd
    does not drop nulls
    :param dataframe: database contaning at least pdb_id, sequences
    :return: new dataframe with the specified characteristics
    """
    if dataframe.empty:
        print('DataFrame is empty!')

    # keep native interactions
    unique_pdb_id = dataframe.filter(items=['Entry','PDB_id','protein_1_seq','protein_2_seq','Wild-type KD (M)'])
    unique_pdb_id.drop_duplicates(inplace=True)
    unique_pdb_id.reset_index(drop=True, inplace=True)
    unique_pdb_id.columns = ['Entry', 'PDB_id', 'protein_1_seq', 'protein_2_seq', 'Kd [M]']
    is_mutated = []
    for j in range(0, len(unique_pdb_id)):
        is_mutated.append('False')
    # make a dataframe for the WT interactions
    wt_data = pd.DataFrame(is_mutated)
    wt_data.columns = ['Is mutant?']
    wt_interactions = pd.concat([unique_pdb_id, wt_data], axis=1)

    # now we make a dataframe with the mutated sequences
    mutated_pdb_id = dataframe.filter(items=['Entry', 'PDB_id', 'mutated_protein_1_seq', 'mutated_protein_2_seq',
                                             'Mutant KD (M)'])
    print(mutated_pdb_id)
    mutated = []
    for i in range(0, len(mutated_pdb_id)):
        mutated.append("True")
    # make a dataframe for the WT interactions
    mut_data = pd.DataFrame(mutated, columns=['Is mutant?'])
    mut_interactions = pd.concat([mutated_pdb_id, mut_data], axis=1)
    print(mut_interactions)
    #mut_interactions.columns = ['Entry', 'PDB_id', 'protein_1_seq', 'protein_2_seq', 'Kd [M]', 'Is mutant?']


    # create global dataframe with both kinds of interactions
    all_interactions = pd.concat([wt_interactions, mut_interactions],axis=0)
    # all_interactions.columns = ['PDB_id', 'protein_1_seq', 'protein_2_seq', 'Kd [M]', 'Is mutant?']
    all_interactions.sort_values(by=['PDB_id'], inplace=True)
    all_interactions.reset_index(drop=True, inplace=True)

    # drop duplicates
    all_interactions.drop_duplicates(subset=['PDB_id', 'protein_1_seq', 'protein_2_seq', 'Is mutant?', 'Kd [M]'], inplace=True)
    all_interactions.reset_index(drop=True, inplace=True)

    return [all_interactions, wt_interactions, mut_interactions]
    # return mut_interactions


def convert(value,unit):
    # takes the value and unit -> value in M
    if unit == 'mM':
        new_value = value * 10e-3
    elif unit == 'uM':
        new_value = value * 10e-6
    elif unit == 'nM':
        new_value = value * 10e-9
    elif unit == 'pM':
        new_value = value * 10e-12
    else:
        new_value = value * 10e-15
    return new_value


def get_unique_seq_pdb(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get the unique sequences from a database DataFrame and add an identifier with numbers for multiple mutant variants.
    :param df: Database DataFrame with columns 'PDB_id', 'protein_1_seq', 'protein_2_seq', 'Is mutant?' at least.
    :return: DataFrame with unique sequences and identifiers.
    """
    df_1 = df[['PDB_id', 'protein_1_seq', 'Is mutant?']].copy()
    df_1.columns = ['PDB_id', 'protein_seq', 'Is mutant?']
    df_1['id'] = df_1['PDB_id'] + '_P1'

    df_2 = df[['PDB_id', 'protein_2_seq', 'Is mutant?']].copy()
    df_2.columns = ['PDB_id', 'protein_seq', 'Is mutant?']
    df_2['id'] = df_2['PDB_id'] + '_P2'

    # Create a defaultdict to store mutant counts for each PDB_id and protein
    mutant_counts = defaultdict(lambda: defaultdict(int))

    # Iterate through the DataFrame to add mutant numbers to identifiers
    for idx, row in df.iterrows():
        pdb_id = row['PDB_id']
        protein = 'P1' if idx % 2 == 0 else 'P2'  # Alternating between P1 and P2
        is_mutant = row['Is mutant?']

        if is_mutant:
            mutant_counts[pdb_id][protein] += 1
            mutant_number = mutant_counts[pdb_id][protein]
            df_1.loc[df_1['PDB_id'] == pdb_id, 'id'] = f'{pdb_id}_{protein}_m{mutant_number}'
            df_2.loc[df_2['PDB_id'] == pdb_id, 'id'] = f'{pdb_id}_{protein}_m{mutant_number}'
        else:
            df_1.loc[df_1['PDB_id'] == pdb_id, 'id'] = f'{pdb_id}_{protein}'
            df_2.loc[df_2['PDB_id'] == pdb_id, 'id'] = f'{pdb_id}_{protein}'

    df_all_seq = pd.concat([df_1, df_2])

    df_unique_seq = df_all_seq.drop_duplicates(subset=['protein_seq', 'Is mutant?', 'id'], ignore_index=True)
    df_uniq = df_unique_seq[['protein_seq', 'id']].copy()
    df_uniq.drop_duplicates(subset=['protein_seq'], inplace=True)
    df_uniq.reset_index(drop=True, inplace=True)

    return df_uniq

