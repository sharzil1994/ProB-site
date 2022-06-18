# ProB-site
# Intro
ProB-site is a CNN model that predicts binding sites of protien protein interactions. It utilizes evolutionalry information and predicted secondary struacture information extracted from protein sequences. 
# System Requirment
This model has been developed in Linux environment with:
* python 3.8.10
* numpy 1.21.5
* pandas 1.4.2
* tensor flow 2.3.0
# Software and Database requirments
To run full version of  the ProB-site, it requires following software to extract features
* [Blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and [UniRef90](https://www.uniprot.org/downloads)
* [HH-suite](https://github.com/soedinglab/hh-suite) and [Uniclust30](https://uniclust.mmseqs.com/)
* [DSSP](https://github.com/cmbi/dssp)
# Built software and Database
*1


# Run ProB-site for prediction
For prediction of binding site in a protein run following command:  
``` python Prediction.py -p 3zeu -c D ```  
Here 3zeu is PDB_ID and D is chain of that PDB_ID, program will download necessary PDB file from online database  
# Data, feature and model
we have provided pre-computed feature and a pretrained model for those interested in reproducing the paper  
List of Dataset used in this research in present in Featrues/data_seq folder
Features are stored in numpy format
Program can be run using pre-computed features without installing Blast+, HH-suite, 
and DSSP software. However binding sites of proteins present in pre-computed features list can be predicted.
