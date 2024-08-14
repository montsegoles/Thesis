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
    df['Has_Noncanonical_Protein_1'] = df['seq_1'].apply(has_noncanonical_db)
    df['Has_Noncanonical_Protein_2'] = df['seq_2'].apply(has_noncanonical_db)
    
    # Create masks for filtering
    mask_protein_1 = ~df['Has_Noncanonical_Protein_1'] & ~df['seq_1'].apply(is_long_sequence_db, args=(max_length,))
    mask_protein_2 = ~df['Has_Noncanonical_Protein_2'] & ~df['seq_2'].apply(is_long_sequence_db, args=(max_length,))
    
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
    df['Has_Noncanonical'] = df['seq'].apply(has_noncanonical)
    canonical_df = df.loc[~df['Has_Noncanonical']]
    filtered_df = canonical_df.loc[~canonical_df['seq'].apply(is_long_sequence, max_length_n=length)]
    filtered_df = filtered_df.drop(columns=['Has_Noncanonical'])  # Drop the 'Has_Noncanonical' column
    return filtered_df