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
* Install DSSP  
  install libraries  
  ```sudo apt-get install libboost-all-dev  
  sudo apt-get install -y libz-dev  
  sudo apt-get install -y libbz2-dev  
  sudo apt-get install -y automake  
  sudo apt-get install -y autotools-dev  
  sudo apt-get install -y autoconf 
  ```
  
  From link [DSSP](https://github.com/cmbi/dssp) download `dssp-3.1.4.tar.gz`  
  unzip it using command  
  ```tar -zxvf dssp-3.1.4.tar.gz```  
  compile the program using follwoing command  
  ```cd dssp-3.1.4  
  ./autogen.sh  
  ./configure  
  make 
  ```
  * Install Blast+ and database
    download [uniref90.fasta.gz](link: https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/)  
    unzip downloaded niref90.fasta.gz file  
    ```gzip -d uniref90.fasta.gz```  
    download [ncbi-blast-2.13.0+-x64-linux.tar.gz](link https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  
    unzip ncbi-blast-2.13.0+-x64-linux.tar.gz  
    ```tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz```  
    In created folder run following command, making sure ./uniref90.fasta proper path is given  
   ```./ncbi-blast-2.13.0+/bin/makeblastdb -in ./uniref90.fasta -parse_seqids -blastdb-version 5 -title "unirefdb" -dbtype prot```
  * Install HH-suite and database
   
   download [uniclust30_2017_10.tar.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/) and 
   [uniclust30_2017_10_hhsuite.tar.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/) and   
   [uniclust_uniprot_mapping.tsv.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/)  

2.unzip downloaded file
  ```tar -zxvf uniref30_2017_10.tar.gz
  tar -zxvf uniclust30_2017_10_hhsuite.tar.gz
  tar -zxvf uniclust_uniprot_mapping.tsv.gz
   ```
   prerequisits, install following dependencies  
  ```sudo apt install pigz
  sudo apt install libopenmpi-dev
  sudo apt install sed
  sudo apt install md5deep
  sudo apt install clustalo
  sudo apt install kalign
  sudo apt install gawk
  sudo apt install node-connect-timeout  
  sudo apt-get install tar 
  ```
 install HH-suite software using follwing commands  
  ```git clone https://github.com/soedinglab/hh-suite.git
  mkdir -p hh-suite/build && cd hh-suite/build
  cmake -DCMAKE_INSTALL_PREFIX=. ..
  make -j 4 && make install
  export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
  ```
  

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
