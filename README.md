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
  ```
  sudo apt-get install libboost-all-dev  
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
  ```
  cd dssp-3.1.4  
  ./autogen.sh  
  ./configure  
  make 
  ```  
  Here dssp-3.1.4/mkdssp is DSSP software path  
* Install Blast+ and database  
  download [uniref90.fasta.gz](link: https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/)  
  unzip downloaded niref90.fasta.gz file  
  ```gzip -d uniref90.fasta.gz```  
  download [ncbi-blast-2.13.0+-x64-linux.tar.gz](link https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  
  unzip ncbi-blast-2.13.0+-x64-linux.tar.gz  
  ```tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz```  
  In created folder run following command, making sure ./uniref90.fasta proper path is given  
  ```./ncbi-blast-2.13.0+/bin/makeblastdb -in ./uniref90.fasta -parse_seqids -blastdb-version 5 -title "unirefdb" -dbtype prot``` 
  The output of makeblastbd command will give the PSIBlast Database  
  Here ncbi-blast-2.10.1+/bin/psiblast is PSIBLAST software path   
  
* Install HH-suite and database  
   
  download [uniclust30_2017_10.tar.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/) and 
  [uniclust30_2017_10_hhsuite.tar.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/) and   
  [uniclust_uniprot_mapping.tsv.gz](link: http://gwdu111.gwdg.de/~compbiol/uniclust/2017_10/)  

   unzip downloaded file
   ```
   tar -zxvf uniref30_2017_10.tar.gz
   tar -zxvf uniclust30_2017_10_hhsuite.tar.gz
   tar -zxvf uniclust_uniprot_mapping.tsv.gz
   ```  
   The output of above commands will give the HH-Suite Database  
   Pre-requisits, install following dependencies  
   ```
   sudo apt install pigz
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
   ```
   git clone https://github.com/soedinglab/hh-suite.git
   mkdir -p hh-suite/build && cd hh-suite/build
   cmake -DCMAKE_INSTALL_PREFIX=. ..
   make -j 4 && make install
   export PATH="$(pwd)/bin:$(pwd)/scripts:$PATH"
   ``` 
   Here hhsuite-3.0.3/build/bin/hhblits is HHsutie software path   
 * Updating software and data bases path in software  
   In software open featuer_extarction.py and update paths according to your device  
   ```
   In line 33 give correct path to data_path='./data_ext/'  
   In line 35 give correct path to dssp = Software_path + "dssp-3.1.4/mkdssp"  
   In line 36 give correct path to PSIBLAST = Software_path + "ncbi-blast-2.13.1+/bin/psiblast"  
   In line 37 give correct path to HHBLITS = Software_path + "hhsuite-3.0.3/build/bin/hhblits"  
   In line 38 give correct path to UR90 = "./unirefdb/uniref90.fasta"   
   In line 39 give correct path to HHDB = "./uniclust30_2017_10"
   ```  
   In software open prediction.py and update paths according to your device  
   ```
   In line 28 give correct path to  pre_path='./Feature/'  
   In line 29 give correct path to fea_path='./data_ext/'  
   ```

# Run ProB-site for prediction
For prediction of binding site in a protein run following command:  
``` python predictor.py -p 3zeu -c D ```  
Here 3zeu is PDB_ID and D is chain of that PDB_ID, program will download necessary PDB file from online database  
# Data, feature and model
we have provided pre-computed feature and a pretrained model for those interested in reproducing the paper  
List of Dataset used in this research in present in Featrues/data_seq folder  
Features are stored in numpy format  
Program can be run using pre-computed features without installing Blast+, HH-suite, 
and DSSP software. However binding sites of proteins present in pre-computed features list can only be predicted
# Cite
For citation use the given BibTeX format  
@article{khan2022prob,  
  title={ProB-Site: Protein Binding Site Prediction Using Local Features},  
  author={Khan, Sharzil Haris and Tayara, Hilal and Chong, Kil To},  
  journal={Cells},  
  volume={11},  
  number={13},  
  pages={2117},  
  year={2022},  
  publisher={MDPI}  
}
