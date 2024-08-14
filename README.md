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
    * [main.py](./Chapter_6/Affinity/main.py): Main script for the processing of SKEMPI, PROXiMATE, and PDBbind datasets.
    * [processing.py](./Chapter_6/Affinity/processing.py): Script used to delete sequences with non-canonical amino acids and those sequences over 1024 amino acids long. 
    * [processing_functions.py](./Chapter_6/Affinity/processing_functions.py): Script containing the developed functions to pre-process the datasets. It includes functions to retrieve protein sequences from the PDB database, to eliminate multimer complexes, to apply mutations, and many more. 
* **Binary**   
    Contains folders which have the processing of the binary datasets. Each subfolder has its specific inputs, outputs, and respective scripts used for processing. These scripts were developed in conjunction with @JulianGarciaVinuesa.
    * [negatome_processing](./Chapter_6/Binary/Binary/Binary_Negatives/negatome_processing.ipynb): Jupyter Notebook with the pricessing procedement for the negatome database negative interaction data.
    * [positives_processing](./Chapter_6/Binary/Binary/Binary_Positives/Positives_Processing.ipynb): Jupyter Notebook with the processing of the positive dataset from Signor, Mint, and Hippie datasets. 
    * [Binary_functions.py](Chapter_6/Binary/Binary/Binary_processing/Binary_functions.py): Script containing the functions to process and create the binary dataset.
    * [Binary_processing](Chapter_6/Binary/Binary/Binary_processing/Binary_processing.ipynb): Jupyter notebook that employs the binary processing functions to obtain the binary interaction dataset. It fuses both positive and negative datasets with an interaction label.
* **Inputs**
    Folder with the input datasets used to build both affinity and binary datasets. If the file was too heavy to be uploaded, a link to download it from Google Drive is available in the [README](./Chapter_6/Inputs/README.txt) file. Source databases:
    * [SKEMPI](https://life.bsc.es/pid/skempi2/database/index)
    * [PROXiMATE](https://www.iitm.ac.in/bioinfo/PROXiMATE/index.html)
    * [PDBbind](http://www.pdbbind.org.cn)
    * [Hippie](https://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)
    * [MINT](https://mint.bio.uniroma2.it)
    * [Signor](https://signor.uniroma2.it)
    * [Negatome](https://mips.helmholtz-muenchen.de/proj/ppi/negatome/)
* **Outputs**
    This folder contains the main outputs of the processing scrips, ie., the preprocessed datasets of binary and affinity interactions. The main files will be briefly described next:
    * [Affinity_database_C.csv](./Chapter_6/Outputs/Affinity_database_C.csv): Dataset of binding affinity expressed in the dissociation constant ([M]). It contains protein pairs up to 1024 amino acid long, and only with canonical amino acids. It was processed from the SKEMPI, PROXiMATE, and PDBbind v2020 datasets. It contains mutated sequences that are labeled with a "1" value in the "Is mutated?" column. The columns of the dataset are: Entry,PDB_id,protein_1_seq,protein_2_seq,Kd [M],Is mutant?,Database. 
    * [Binary_database_C.csv](https://drive.google.com/file/d/1Dwe5p1f6wILITPlF3Qe5KhrTcn-l2TD_/view?usp=drive_link): Dataset of the binary interactions, the interacting or non-interacting protein pairs. It contains protein pairs up to 1024 amino acid long, and only with canonical amino acids. It was created from the Hippie, MINT, Signor, and Negatome datasets. 
* **dataset_construction**
    Contains a script to create the datasets to train a machine learning or deep learning model, from the binary and affinity datasets. It separates each input dataset onto a training, testing, validation, and benchmark dataset. This folder also features a Jupyter Notebook with the analysis of some of the datasets available in the [Outputs](./Chapter_6/Outputs/) folder.

Contact: montserrat.goles@ug.uchile.cl