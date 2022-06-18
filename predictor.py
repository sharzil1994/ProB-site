# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 06:36:25 2022

@author: PC
"""

import pickle
import numpy as np

# from tqdm import tqdm
# from Bio import SeqIO
import os
import pandas as pd

import subprocess
from feature_extarction import *
import argparse


# os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID";
# # The GPU id to use, usually either "0" or "1";
# os.environ["CUDA_VISIBLE_DEVICES"] = "1";


def site_predictor(input_ID):
    error_code_dic = {1:"PDB not exist", 2:"chain not exist", 3: "DSSP too long",4: "Fail to pad DSSP",5:"software path error", 6:"database path error" }
    pre_path='./Feature/' # change path according to your device(pre computed features)
    fea_path='./data_ext/' #change path according to your device(path of computed features)
    seq=""
    path=pre_path
    with open(f'{path}/data_seq/data_ID.dat','rb') as FH:
        data0 =pickle.load(FH)
        
        
    pssm=np.load(f'{path}pssm/pssm.npy',allow_pickle=True)
    dssp=np.load(f'{path}dssp/dssp.npy',allow_pickle=True)
    hhm=np.load(f'{path}hhm/hhm.npy',allow_pickle=True)
    
    found=0
    i_0=0
    for index_0 in data0[:,0]:  
        if input_ID.lower()==index_0.lower():
            found=1
            seq=data0[i_0,1]
            tEMP_0=pssm[1,i_0]
            tEMP_1=dssp[1,i_0]
            tEMP_2=hhm[1,i_0]
            break
        i_0=i_0+1

    
    if found==1:
        pass
    else :
        print("Sequence feature not present in data base, extracting features, it will take several minutes")
        error,seq=feature_extraction(input_ID)
        
        if error!=0:
            return error_code_dic[error],'',''
        
        path=fea_path
        tEMP_0=np.load(f'{path}pssm/{input_ID}.npy')
        tEMP_1=np.load(f'{path}dssp/{input_ID}.npy')
        tEMP_2=np.load(f'{path}hhm/{input_ID}.npy')
    
    
    def windowing(features,seq,w_size,start,stop):   
    
        # all_features=features
        seq_len=len(seq)
        
        
        # a=all_features
        fea_len=np.shape(features)[1]
        # print(fea_len)
    
        finalout1=np.zeros([seq_len,w_size,fea_len],'float')
        # finalout1=np.zeros([to,w_size,69],'float')
        
        l=0
        for j in range( 0, seq_len):
            
            for k in range(0,w_size):
                
                k1=int(j+k-((w_size-1)/2))
                
                if k1<0 or k1 > seq_len-1:
                    pass
                else:
                    finalout1[l,k,:]=features[k1,start:stop]
    
            l=l+1    
        return finalout1
    
    import tensorflow as tf
    
    w_size=19
    f1=windowing(tEMP_0,seq,w_size,0,20)
    f2=windowing(tEMP_1,seq,w_size,0,14)
    f3=windowing(tEMP_2,seq,w_size,0,20)
    new_model = tf.keras.models.load_model('./model/model_19.h5')
    pred_y = new_model.predict([f1,f2,f3])
    
    for i in range(len(pred_y)):
        if pred_y[i] < 0.5:
            pred_y[i] = 0;
        else:
            pred_y[i] = 1;
    
    pred_y=np.squeeze(pred_y)        
    pred_y=pred_y.astype('int')
    
    out=''
    for index_0 in pred_y:
        out=out+str(index_0)        

    return "", seq , out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--protein", type = str, help = "PDBID (e.g. 3zeu)")
    parser.add_argument("-c", "--chain", type = str, help = "chain identifier(e.g D)")
    args = parser.parse_args()
    if args.chain == None:
        print("Chain identifier is not provided!")
    elif args.chain not in list(string.ascii_letters + string.digits):
        print("Invalid chain identifier!")
    else: # input by PDBID
        if args.protein == None or len(args.protein) != 4:
            print("Invalid PDB ID!")
        else:
            print(args.protein + args.chain)
            error,seq,pred=site_predictor(args.protein + args.chain )
    if error !='':
        print(error)
    else:
        print(seq)
        print(pred)
