import pandas as pd
import requests
import numpy as np
from processing_functions import *

# This script is to get the unique sequences of the complexes from the databeses
# we import the databases
df = pd.read_csv('database_processing\Outputs\Affinity_database.csv')
# %%
unique_seq = get_unique_seq_pdb(df)

# %% Second, we cut those with non-canonical amino acids and length over the maximum
df_cut = filter_canonical_and_length(1024, unique_seq)
df_cut.to_csv(r'database_processing\Outputs\unique_affinity_seq.csv', index=False)
#%% Now we get a sample

sample_aff = df_cut.sample(n=100)
#%%
sample_aff.to_csv(r'database_processing\Outputs\Affinity_sample.csv', index=False)
#%% 
df_affinity_canonical = filter_canonical_and_length_db(1024,df)
df_affinity_canonical.to_csv(r'database_processing\Outputs\Affinity_database_C.csv', index=False)