
from typing import List, Any

import pandas as pd
import numpy as np

from processing_functions import *
#%%
# First Proximate
df_proximate = pd.read_csv("Inputs/proximate-20mar2021.csv", encoding='windows-1252')
df_proximate.columns = ['Entry', 'PDB_id', 'Interacting Chains', 'Mutation(s)',
                        'Protein 1', 'Protein 2', 'Functional Class', 'Technique used', 'Ionic Conditions',
                        'Temperature (K)', 'pH', 'Database Source', 'PMID/Reference', 'Wild-type KD (M)',
                        'Mutant KD (M)',
                        'Wild-type kon (1/Ms)', 'Mutant kon (1/Ms)', 'Wild-type koff (1/s)', 'Mutant koff (1/s)',
                        'Wild-type ΔG (kcal/mol)', 'ΔΔG (kcal/mol)', 'Wild-type ΔH (kcal/mol)', 'Mutant ΔH (kcal/mol)',
                        'Wild-type ΔS (cal/mol K)', 'Mutant ΔS (cal/mol K)', 'Notes']
#%%
# separate chains
df_proximate[['Chain 1', 'Chain 2']] = df_proximate['Interacting Chains'].str.split(':', expand=True)
# drop empty chains
proximate = drop_nulls_by('Interacting Chains', df_proximate)
proximate.reset_index(drop=True, inplace=True)
proximate_dimers = proximate.copy()
# drop nulls by
proximate_no_nulls = drop_nulls_by('PDB_id', proximate_dimers)
proximate_no_nulls2 = drop_nulls_by('Wild-type KD (M)',proximate_no_nulls)
#%%
# process daataframe: eliminate multimers and get sequences
proximate_dimers_seq = process_dataframe(proximate_no_nulls2)


#%%

proximate_no_mut = proximate_dimers_seq.copy()
#%%
proximate_mut = apply_mutations(proximate_no_mut,'Chain 1','Chain 2', 'Mutation(s)')
#%%
proximate_processed = process_mutated_db(proximate_mut)

#%%
proximate_WT = proximate_processed[1].copy()
proximate_mut = proximate_processed[2].copy()

#%%
proximate_WT.drop_duplicates(subset=['PDB_id','Kd [M]'],inplace=True)
#%%
proximate_mut.drop_duplicates(subset=['PDB_id','mutated_protein_1_seq','mutated_protein_2_seq'],inplace=True)
proximate_mut.dropna(subset=['Mutant KD (M)'],inplace=True)
#%%
proximate_WT.reset_index(drop=True,inplace=True)
proximate_mut.reset_index(drop=True,inplace=True)
#%%

proximate_mut.columns = ['Entry', 'PDB_id', 'protein_1_seq', 'protein_2_seq', 'Kd [M]', 'Is mutant?']

#%% merge both dataframes
proximate_both = pd.concat([proximate_WT, proximate_mut], axis=0)
#%%
proximate_both['Database'] = 'PROXIMATE'

# proximate_both.to_csv('Outputs/PROXIMATE_DIMERS.csv',index=False)
#%%
PROXIMATE_ND = drop_duplicates_by_sequence(proximate_both)

PROXIMATE_ND.to_csv('Outputs/PROXIMATE_VPDB.csv',index=False)

#%%
# for skempi

all_skempi = pd.read_csv("Inputs/skempi_v2.csv", sep=';')
all_skempi.columns = ['PDB_id', 'Mutation(s)_PDB', 'Mutation(s)_cleaned', 'iMutation_Location(s)', 'Hold_out_type',
                      'Hold_out_proteins', 'Affinity_mut (M)', 'Affinity_mut_parsed', 'Wild-type KD (M)',
                      'Affinity_wt_parsed', 'Reference', 'Protein 1', 'Protein 2', 'Temperature',
                      'kon_mut (M^(-1)s^(-1))', 'kon_mut_parsed', 'kon_wt (M^(-1)s^(-1))', 'kon_wt_parsed',
                      'koff_mut (s^(-1))', 'koff_mut_parsed', 'koff_wt (s^(-1))', 'koff_wt_parsed',
                      'dH_mut (kcal mol^(-1))', 'dH_wt (kcal mol^(-1))', 'dS_mut (cal mol^(-1) K^(-1))',
                      'dS_wt (cal mol^(-1) K^(-1))', 'Notes', 'Method', 'SKEMPI version']
all_skempi[['PDB_id', 'Chain 1', 'Chain 2']] = all_skempi['PDB_id'].str.split('_', expand=True)
all_skempi['Entry'] = range(1, len(all_skempi) + 1)
skempi_all_len = len(all_skempi)
skempi = drop_nulls_by('PDB_id', all_skempi)
skempi.reset_index(drop=True, inplace=True)
# drop nulls by chains
skempi_with_nulls = drop_nulls_by('Chain 1', skempi)
skempi_with_nulls.reset_index(drop=True, inplace=True)
skempi_no_nulls = drop_nulls_by('Chain 2', skempi)
skempi_no_nulls.reset_index(drop=True, inplace=True)
#%%
# eliminate multimers and get sequences

skempi_seq = process_dataframe(skempi_no_nulls)

#%%
skempi_no_mut = skempi_seq.copy()
#%%
skempi_no_mut.reset_index(drop=True, inplace=True)
#%%
skempi_no_mut.drop([929],inplace=True)
skempi_no_mut.drop([930],inplace=True)
skempi_no_mut.drop([931],inplace=True)
skempi_no_mut.drop([932],inplace=True)
skempi_no_mut.drop([933],inplace=True)
skempi_no_mut.drop([934],inplace=True)
skempi_no_mut.drop([950],inplace=True)
skempi_no_mut.reset_index(drop=True, inplace= True)
#%%
skempi_no_mut['Mutant KD (M)'] = skempi_no_mut['Affinity_mut (M)']
#%%
skempi_mut = apply_mutations(skempi_no_mut,'Chain 1','Chain 2', 'Mutation(s)_PDB')

#%%
skempi_processed = process_mutated_db(skempi_mut)

#%%
skempi_WT = skempi_processed[1].copy()
skempi_mut = skempi_processed[2].copy()

#%%
skempi_WT.drop_duplicates(subset=['PDB_id', 'Kd [M]'],inplace=True)
#%%
skempi_mut.drop_duplicates(subset=['PDB_id', 'mutated_protein_1_seq', 'mutated_protein_2_seq'],inplace=True)
skempi_mut.dropna(subset=['Mutant KD (M)'], inplace=True)
#%%
skempi_WT.reset_index(drop=True, inplace=True)
skempi_mut.reset_index(drop=True, inplace=True)
#%%

skempi_mut.columns = ['Entry', 'PDB_id', 'protein_1_seq', 'protein_2_seq', 'Kd [M]', 'Is mutant?']

#%% merge both dataframes
skempi_both = pd.concat([skempi_WT, skempi_mut], axis=0)
#%%
skempi_both['Database'] = 'SKEMPI'

# skempi_both.to_csv('Outputs/SKEMPI_DIMERS.csv', index=False)


#%%
skempi_both.info()
SKEMPI_ND = drop_duplicates_by_sequence(skempi_both)

#%%

SKEMPI_ND.to_csv('Outputs/SKEMPI_VPDB.csv',index=False)
#%% Now PDBbind
all_pdbbind = pd.read_table("Inputs/INDEX_general_PP_text.txt")
all_pdbbind.columns = ['all']
all_pdbbind[['data','reference and notes']] = all_pdbbind['all'].str.split('//',expand=True)
all_pdbbind['Entry'] = range(1, len(all_pdbbind) + 1)
all_pdbbind.drop(columns=['all'], inplace=True)
all_pdbbind[['PDB_id','resolution','release year','binding data','none1','none2','none3','none4']] =all_pdbbind['data'].str.split('  ',expand=True)
all_pdbbind.drop(columns=['none1','none2','none3','none4','data'],inplace=True)
pdb_bind = all_pdbbind.filter(items=['Entry','PDB_id','binding data'])
# separate constant type from values
pdb_bind[['constant','value']] = pdb_bind['binding data'].str.split('=',expand=True)
#%% filter non-Kd constants
for i in range(0,len(pdb_bind)):
    if pdb_bind['constant'][i] != 'Kd':
        pdb_bind.drop([i],inplace=True)
pdb_bind.reset_index(inplace=True)
#%% separate value column
db_split = pdb_bind.filter(items=['Entry','PDB_id'])
db_split['value'] = ''
d_value = []
db_split['unit'] = ''
d_unit = []
for i in range(0,len(pdb_bind)):
    estring = pdb_bind['value'][i]
    length = len(estring)
    # get the value
    d_value.append(estring[0:length-2])
    # get the dimention
    d_unit.append(estring[-2:])
db_split['value'] = d_value
db_split['unit'] = d_unit
db_split['value'] = pd.to_numeric(db_split['value'])
#%%
pdbbind_new = db_split.filter(items=['Entry','PDB_id'])
pdbbind_db = pdbbind_new.copy()
pdbbind_db['Kd [M]'] = ''
lista_values = []
for i in range(0,len(pdb_bind)):
    value = db_split['value'][i]
    unit = db_split['unit'][i]
    lista_values.append(convert(value, unit))
pdbbind_db['Kd [M]'] = lista_values
pdbbind_db['Kd [M]'] = pd.to_numeric(pdbbind_db['Kd [M]'])
#%% get sequences and drop multimers
df_pdbbind = process_dataframe(pdbbind_db)
df_pdbbind.reset_index(drop=True,inplace=True)
#%%
df_pdbbind.to_csv('PDBbind_vPDB.csv',index=False)
#%% Now we merge

proximate = pd.read_csv('Outputs/PROXIMATE_VPDB.csv')
skempi = pd.read_csv('Outputs/SKEMPI_VPDB.csv')
PDBbind = df_pdbbind.copy()

# add the database column to PDBbind
PDBbind['Database'] = 'PDBbind'
PDBbind['Is mutant?'] = 'False'

#%% drop sorted columns
proximate.drop(columns=['protein_1_seq_sorted','protein_2_seq_sorted'], inplace=True)
#%%
skempi.drop(columns=['protein_1_seq_sorted','protein_2_seq_sorted'],inplace=True)
#%% now me concat
df1 = pd.concat([proximate,skempi],axis=0)
#%%
affinity_database_M = pd.concat([df1, PDBbind],axis=0)
#%%
duplicados = affinity_database_M[affinity_database_M.duplicated(subset=['PDB_id','protein_1_seq','protein_2_seq'])]
if duplicados.empty:
    print("No duplicate interactions found.")
else:
    print("Duplicate interactions found:")
    print(duplicados)
#%% sin duplicados

affinity_database = affinity_database_M.drop_duplicates(subset=['PDB_id','protein_1_seq','protein_2_seq','Is mutant?'])
#%% export

affinity_database.to_csv('Outputs/Affinity_mutants.csv',index=False)
#%%
affinity_database.reset_index(drop=True, inplace=True)
# check if Kd values are valid
# affinity_database['Kd [M]'] = pd.to_numeric(affinity_database['Kd [M]'])
#%%
without_nb = affinity_database.copy()
for i in range(0, len(affinity_database)):
    if affinity_database['Kd [M]'][i] == 'n.b':
        print(affinity_database['Kd [M]'][i])
        without_nb.drop([i], inplace=True)
without_nb.reset_index(drop=True,inplace=True)
#%% to numeric
# without_nb['Kd [M]'] = pd.to_numeric(without_nb['Kd [M]'])
without_extras = without_nb.copy()
#%%
for i in range(0, len(without_nb)):
    if without_nb['Kd [M]'][i] == "unf":
        print(without_nb['Kd [M]'][i])
        without_extras.drop([i], inplace=True)
    elif '~' in str(without_nb['Kd [M]'][i]):
        print('str ~', without_nb['Kd [M]'][i])
        without_extras.drop([i], inplace=True)
    elif '>' in str(without_nb['Kd [M]'][i]):
        print('str >', without_nb['Kd [M]'][i])
        without_extras.drop([i], inplace=True)
    elif '<' in str(without_nb['Kd [M]'][i]):
        print('str <', without_nb['Kd [M]'][i])
        without_extras.drop([i], inplace=True)
    elif 'n.b.' in str(without_nb['Kd [M]'][i]):
        print('str n.b.', without_nb['Kd [M]'][i])
        without_extras.drop([i], inplace=True)
without_extras.reset_index(drop=True,inplace=True)
#%%
without_extras['Kd [M]'] = pd.to_numeric(without_extras['Kd [M]'])
#%% drop duplicate interactions by protein sequences
final_db = drop_duplicates_by_sequence(without_extras)
final_db.drop(columns=['protein_1_seq_sorted', 'protein_2_seq_sorted'],inplace=True)
final_db.to_csv('Outputs/Affinity_database.csv', index=False)