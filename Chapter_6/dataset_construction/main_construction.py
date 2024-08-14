#%%
# import libraries
# pip install sklearn 
import pandas as pd
import sklearn as sk
from sklearn.model_selection import train_test_split
# import datasets

# Affinnity 
affinity_all = pd.read_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\Outputs\Affinity_database_C.csv')
#affinity_all.columns = ['Entry','PDB_id','protein_1_seq','protein_2_seq','Kd [M]','Is mutant?','Database']

# separate by mutant or not
mutant_df = affinity_all[affinity_all['Is mutant?'] == True]
non_mutant_df = affinity_all[affinity_all['Is mutant?'] == False]

#%%
# Binary
binary_df = pd.read_csv(r"C:\Users\monts\Downloads\Binary_database_C.csv")
#binary_df.columns = ['id_1','seq_1','id_2','seq_2','dataset','interaction']
binary_df.rename(columns={"Interaccion": "interaction"}, inplace=True)

#%%
# separate binary by interacting or not interacting
binary_int = binary_df[binary_df['interaction'] == 1]
binary_non = binary_df[binary_df['interaction'] == 0]

#%% create benchmark dataset 
'''
Benchmark dataset
10% of each:
    - 10% positive = 627479*10% =  62747.9  approx 62,748 # len(binary_int)
    - 10% negative = 1965*10% =  196.5 approx 197 # len(binary_non)

the divided datasets
    - positive 62,748 / 564731
    - negative 197 / 1768

'''
# Binary

bin_pos_rest, bin_pos_bench = train_test_split(binary_int, test_size=0.1, train_size= 0.9, random_state= 42)

bin_neg_rest, bin_neg_bench = train_test_split(binary_non, test_size=0.1, train_size=0.9, random_state=42)

#%% Affinity (does not have negatives)


# All affinity
aff_rest, aff_bench = train_test_split(affinity_all, test_size=0.1, train_size=0.9, random_state=42)

# mutated only
mut_rest, mut_bench = train_test_split(mutant_df, test_size=0.1, train_size=0.9, random_state=42)

# wild type only (non mutated)
no_mut_rest, no_mut_bench = train_test_split(non_mutant_df,test_size=0.1, train_size=0.9, random_state=42)


#%%
'''
Training set
70% of what was not used for the benchmark dataset

'''
# Binary
bin_pos_train, bin_pos_no_train = train_test_split(bin_pos_rest, test_size=0.3, train_size= 0.7, random_state= 42)

bin_neg_train, bin_neg_no_train  = train_test_split(bin_neg_rest, test_size=0.3, train_size= 0.7, random_state= 42)

# Affinity
aff_train, aff_no_train = train_test_split(aff_rest, test_size=0.3, train_size= 0.7, random_state= 42)

# mutated only
mut_train, mut_no_train = train_test_split(mut_rest, test_size=0.3, train_size= 0.7, random_state= 42)

# wild type only (non mutated)
no_mut_train, no_mut_no_train = train_test_split(no_mut_rest, test_size=0.3, train_size= 0.7, random_state= 42)

#%% 
'''
Validation set
20% of the rest
and 
Testing set 
10% of the rest

'''
# Binary
bin_pos_val, bin_pos_test, = train_test_split(bin_pos_no_train, test_size=1/3, train_size= 2/3, random_state= 42)

bin_neg_val, bin_neg_test =  train_test_split(bin_neg_no_train, test_size=1/3, train_size= 2/3, random_state= 42)

# Affinity
aff_val, aff_test = train_test_split(aff_no_train, test_size=1/3, train_size= 2/3, random_state= 42)

# mutated only 
mut_val, mut_test = train_test_split(mut_no_train, test_size=1/3, train_size= 2/3, random_state= 42)

# wild type only (non mutated)
no_mut_val, no_mut_test = train_test_split(no_mut_no_train, test_size=1/3, train_size= 2/3, random_state= 42)

#%% Save Datasets
'''
benchmark
train
validation
test
for: positive, negative, all affinity, non mutated and mutated

'''
#%%+ Save each dataset as CSV

# Benchmark datasets
bin_pos_bench.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_pos_bench.csv', index=False)
bin_neg_bench.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_neg_bench.csv', index=False)
aff_bench.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\aff_bench.csv', index=False)
mut_bench.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\mut_bench.csv', index=False)
no_mut_bench.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\no_mut_bench.csv', index=False)

# Training datasets
bin_pos_train.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_pos_train.csv', index=False)
bin_neg_train.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_neg_train.csv', index=False)
aff_train.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\aff_train.csv', index=False)
mut_train.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\mut_train.csv', index=False)
no_mut_train.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\no_mut_train.csv', index=False)

# Validation datasets
bin_pos_val.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_pos_val.csv', index=False)
bin_neg_val.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_neg_val.csv', index=False)
aff_val.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\aff_val.csv', index=False)
mut_val.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\mut_val.csv', index=False)
no_mut_val.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\no_mut_val.csv', index=False)

# Testing datasets
bin_pos_test.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_pos_test.csv', index=False)
bin_neg_test.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\bin_neg_test.csv', index=False)
aff_test.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\aff_test.csv', index=False)
mut_test.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\mut_test.csv', index=False)
no_mut_test.to_csv(r'C:\Users\monts\Documents\GitHub\PPI-Project\database_processing\dataset_construction\Outputs\no_mut_test.csv', index=False)
