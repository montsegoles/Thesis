# Thesis
Code and datasets for my thesis: De novo therapeutic peptide design using deep generative models.
---

Chapter 6 contains:
* Code for the processing of affinity and binary PPI datasets.
* Input and output datasets. If the datasets were too heavy to upload to this repository, a README file in each section has links to a Google Drive folder with the remaining datasets.
* Code for the generation of training, testing, validation, and benchmark datasets for developing a PPI DL-based model.
* A text file with additional information about the structure and data of each folder. 

---
## Contents of Chapter 6 Folder


* **Affinity**
    Contains the scripts used to process the binding affinity datasets. It has three files:
    * [main](/Chapter 6/Affinity/main.py): Main script for the processing of SKEMPI, PROXiMATE, and PDBbind datasets.
    * [processing](/Chapter 6/Affinity/processing.py): Script used to delete sequences with non-canonical amino acids and those sequences over 1024 amino acids long. 
    * [processing_functions.py](/Chapter 6/Affinity/processing_functions.py): Script containing the developed functions to pre-process the datasets. It includes functions to retrieve protein sequences from the PDB database, to eliminate multimer complexes, to apply mutations, and many more. 
