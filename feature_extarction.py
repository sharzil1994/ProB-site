#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This code has been developed using Yuan et al. codes
# link of the source code is "https://github.com/biomed-AI/GraphPPIS/blob/232b4b1cfd8b41607c80f4c0b2a1a8627d0bdab8/GraphPPIS_predict.py#L34"

"""
Created on Mon Jun 13 15:02:05 2022

@author: haris
"""

import os, pickle, datetime, argparse, string
import numpy as np
import pandas as pd



def feature_extraction(Input_ID):  
    AA_SYM = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
          "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
    AA_AABR = [index for index in "ACDEFGHIKLMNPQRSTVWY"]
    aa_dict = dict(zip(AA_SYM, AA_AABR))
    
    Max_pssm = np.array([8, 9, 9, 9, 12, 9, 8, 8, 12, 9, 7, 9, 11, 10, 9, 8, 8, 13, 10, 8])
    Min_pssm = np.array([-11,-12,-13,-13,-12,-13,-13,-12,-13,-13,-13,-13,-12,-12,-13,-12,-12,-13,-13,-12])
    Max_hhm = np.array([10655,12141,12162,11354,11802,11835,11457,11686,11806,11262,11571,11979,12234,11884,11732,11508,11207,11388,12201,11743])
    Min_hhm = np.zeros(20)
    
    
    error_code_dic = {"PDB not exist": 1, "chain not exist": 2,  "DSSP too long": 3, "Fail to pad DSSP": 4,"software path error":5, "database path error" : 6}
    
    #Add path here accoding to your device
    data_path='./data_ext/' #output folder path for features
    Software_path = "./"    
    dssp = Software_path + "dssp-3.1.4/mkdssp" #dssp software folder path for features extraction
    PSIBLAST = Software_path + "ncbi-blast-2.13.1+/bin/psiblast" #PSIBLAST software folder path for features extraction
    HHBLITS = Software_path + "hhsuite-3.0.3/build/bin/hhblits" #HH-SUITE software folder path for features extraction
    UR90 = "./unirefdb/uniref90.fasta" #database for pssm path #PSIBLAST DATABASE folder path for features extraction
    HHDB = "./uniclust30_2017_10" #database hmm path #HH-SUITE DATABASE folder path for features extraction
    
    if os.path.exists(dssp) == False:
        return  error_code_dic["software path error"],""

    if os.path.exists(PSIBLAST) == False:
        return  error_code_dic["software path error"],""    
    
    if os.path.exists(HHBLITS) == False:
        return  error_code_dic["software path error"],""
    
    if os.path.exists(UR90) == False:
        return  error_code_dic["database path error"],""
   
    if os.path.exists(HHDB) == False:
        return  error_code_dic["database path error"] ,""  
    
    def pdb_extract(datapath,pdbid,chain):
        ID=pdbid + chain
        if os.path.exists(data_path + "{}.pdb.gz".format(pdbid)) == False:
            os.system("wget -P {} http://www.rcsb.org/pdb/files/{}.pdb.gz".format(data_path, pdbid))
        
        
        if os.path.exists(data_path + "{}.pdb.gz".format(pdbid)) == False:
            return "", error_code_dic["PDB not exist"]
        
        os.system("perl getchain.pl {} {}".format(data_path, ID))
        
        os.system("mv {} {}".format(ID, data_path))
        seq_p = ""
        indicator = -1000
        with open(data_path + ID, "r") as f:
            out = f.readlines()
        for index in out:
            if index[0:4].strip() == "ATOM" and int(index[22:26].strip()) != indicator:
                aa_type = index[17:20].strip()
                seq_p += aa_dict[aa_type]
                indicator=  int(index[22:26].strip())
          
        with open(data_path + ID + ".fa", "w") as f:
            f.write(">" + ID + "\n" + seq_p)
    
        
        if seq_p == "":
            return "", error_code_dic["chain not exist"]
        else:
            return seq_p, 0
    
    
         
    
    def process_dssp(dssp_file):
            
        aa_type = "ACDEFGHIKLMNPQRSTVWY"
        SS_type = "HBEGITSC"
        rASA_std = [115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
                    185, 160, 145, 180, 225, 115, 140, 155, 255, 230]
        
        
        
        
        with open(dssp_file, "r") as f:
            out = f.readlines()
        
        seq_d = ""
        dssp_feature = []
        
        
        p = 0
        while out[p].strip()[0] != "#":
            p += 1
        for i in range(p + 1, len(out)):
            aa = out[i][13]
            if aa == "!" or aa == "*":
                continue
            seq_d += aa
            SS = out[i][16]
            if SS == " ":
                SS = "C"
            SS_vec = np.zeros(9) # The last dim represents "Unknown" for missing residues
            SS_vec[SS_type.find(SS)] = 1
            PHI = float(out[i][103:109].strip())
            PSI = float(out[i][109:115].strip())
            ACC = float(out[i][34:38].strip())
            ASA = min(100, round(ACC / rASA_std[aa_type.find(aa)] * 100)) / 100
            dssp_feature.append(np.concatenate((np.array([PHI, PSI, ASA]), SS_vec)))
        
        dssp_f=np.array(dssp_feature)
        return seq_d,dssp_f
    
    
    
    
    
    
    def pad_dssp(seq, feature, ref_seq): # ref_seq is longer
        padded_feature = []
        SS_vec = np.zeros(9) # The last dim represent "Unknown" for missing residues
        SS_vec[-1] = 1
        padded_item = np.concatenate((np.array([360, 360, 0]), SS_vec))
    
        p_ref = 0
        for i in range(len(seq)):
            while p_ref < len(ref_seq) and seq[i] != ref_seq[p_ref]:
                padded_feature.append(padded_item)
                p_ref += 1
            if p_ref < len(ref_seq): # aa matched
                padded_feature.append(feature[i])
                p_ref += 1
            else: # miss match!
                return np.array([])
    
        if len(padded_feature) != len(ref_seq):
            for i in range(len(ref_seq) - len(padded_feature)):
                padded_feature.append(padded_item)
    
        return np.array(padded_feature)
    
    
    def transform_dssp(dssp_feature):
        angle = dssp_feature[:,0:2]
        ASA_SS = dssp_feature[:,2:]
    
        radian = angle * (np.pi / 180)
        dssp_feature = np.concatenate([np.sin(radian), np.cos(radian), ASA_SS], axis = 1)
    
        return dssp_feature
    
    def get_dssp(ID, PDB_seq, data_path):
        if os.path.exists(data_path + "dssp/{}.dssp".format(ID)) == False:
            os.system("{} -i {} -o {}.dssp".format(dssp, data_path + ID, data_path + "dssp/"  + ID))
        dssp_seq, dssp_matrix = process_dssp(data_path+"dssp/" +  ID + ".dssp")
        if len(dssp_seq) > len(PDB_seq):
            return error_code_dic["DSSP too long"]
        elif len(dssp_seq) < len(PDB_seq):
            padded_dssp_matrix = pad_dssp(dssp_seq, dssp_matrix, PDB_seq)
            if len(padded_dssp_matrix) == 0:
                return error_code_dic["Fail to pad DSSP"]
            else:
                np.save(data_path + "dssp/" + ID, transform_dssp(padded_dssp_matrix))
        else:
            np.save(data_path+ "dssp/" + ID, transform_dssp(dssp_matrix))
        return 0
    
    

    
    
    
    
    def process_pssm(pssm_file):
        with open(pssm_file, "r") as f:
            lines = f.readlines()
        pssm_feature = []
        for line in lines:
            if line == "\n":
                continue
            record = line.strip().split()
            if record[0].isdigit():
                pssm_feature.append([int(x) for x in record[2:22]])
        pssm_feature = (np.array(pssm_feature) - Min_pssm) / (Max_pssm - Min_pssm)
    
        return pssm_feature
    
    
    def process_hhm(hhm_file):
        with open(hhm_file, "r") as f:
            lines = f.readlines()
        hhm_feature = []
        p = 0
        while lines[p][0] != "#":
            p += 1
        p += 5
        for i in range(p, len(lines), 3):
            if lines[i] == "//\n":
                continue
            feature = []
            record = lines[i].strip().split()[2:-1]
            for x in record:
                if x == "*":
                    feature.append(9999)
                else:
                    feature.append(int(x))
            hhm_feature.append(feature)
        hhm_feature = (np.array(hhm_feature) - Min_hhm) / (Max_hhm - Min_hhm)
    
        return hhm_feature
    

    pdbid=Input_ID[0:len(Input_ID)-1]
    chain=Input_ID[len(Input_ID)-1]


    
    
    seq_p,error=pdb_extract(data_path  ,pdbid , chain)
    if error!=0:
        return error,""
    
    
    error=get_dssp(Input_ID, seq_p, data_path)
    
    if error!=0:
        return error,""
    if os.path.exists(data_path + "pssm/{}.pssm".format(Input_ID)) == False:
        os.system("{0} -db {1} -num_iterations 3 -num_alignments 1 -num_threads 2 -query {3}{2}.fa -out {3}{2}.bla -out_ascii_pssm {3}pssm/{2}.pssm".format(PSIBLAST, UR90, Input_ID, data_path))
    if os.path.exists(data_path + "hhm/{}.hhm".format(Input_ID)) == False:
        os.system("{0} -i {2}{1}.fa -ohhm {2}hhm/{1}.hhm -oa3m {2}{1}.a3m -d {3} -v 0 -maxres 40000 -cpu 6 -Z 0 -o {2}{1}.hhr".format(HHBLITS, Input_ID, data_path, HHDB))
    
    pssm_matrix = process_pssm(data_path + "pssm/" + Input_ID + ".pssm")
    
    hhm_matrix = process_hhm(data_path + "hhm/" + Input_ID + ".hhm")
    
    
    
    np.save(data_path + "pssm/" + Input_ID, pssm_matrix)
    np.save(data_path + "hhm/" + Input_ID, hhm_matrix)
    
    return error,seq_p





