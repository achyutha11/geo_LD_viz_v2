#!python 

"""
    Functions to convert K-banded and adaptive-K matrices to P x P matrix (P is number of SNPs)
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def kxp_convert(kxpmat):
    """
    Convert K x P matrix to P x P matrix where P is the number of SNPs
    Argument:
    kxpmat : K x P matrix 
    """

    K,P = kxpmat.shape
    assert(K > 1)
    assert(P > K)
    
    # Create zeros matrix (nsnp x nsnp)
    emp = np.zeros((P,P))

    # Iterate through columns of K x P matrix, appending each to appropriate column in zeros matrix
    for i in tqdm(range(P)):
        if i < (P - K):     
            emp[i+1:i+1+kxpmat[2:,i].shape[0],i] = kxpmat[2:,i]
        # Accounting for edge cases
        else:
            emp[i+1:,i] = kxpmat[1:,i][:kxpmat.shape[1]-i-1]
    
    return emp

def kxp_file_convert(kxp_file):
    """
    Convert K x P file to P x P matrix where P is the number of SNPs
    Argument:
    kxp_file : K x P npz file path
    """
  
    inter = np.load(kxp_file)
    kxpmat = inter['ld_mat']
    
    return kxp_convert(kxpmat)

def ragged_convert(array_list,blen):
    """
    Convert list of adaptive-K arrays to P x P matrix where P is the number of SNPs
    Arguments:
    array_list : List of adaptive-K arrays
    blen : Block length used in the function to generate adaptive arrays
    """
   
    assert(blen < len(array_list))
    assert(len(array_list) > 1)
    P = len(array_list) 
    
    # Creating empty matrix
    mat = np.zeros((P+1,P+1))
    
    # Iterating list of arrays and appending each list to appropriate position of P x P matrix
    for i in tqdm(range(P)):
        if i < P - blen: 
            if array_list[i].shape[0] >= blen:
                mat[i+1:array_list[i].shape[0]+i+1,i] = array_list[i]
            elif array_list[i].shape[0] == 0:
                continue
            elif array_list[i].shape[0] > 0 and array_list[i].shape[0] < blen:
                if all(array_list[i][m].shape[0]%blen==0 for m in list(range(array_list[i].shape[0]))):
                    index = array_list[i].shape[0]*blen
                    mat[i+1:index+i+1,i] = np.concatenate(array_list[i]) 
                else:
                    lengths = []
                    for k in range(array_list[i].shape[0]):
                        lengths.append(array_list[i][k].shape[0])
                    index = sum(lengths)
                    mat[i+1:index+i+1,i] = np.concatenate(array_list[i]) 
        else:
            mat[i+1:array_list[i].shape[0]+i+1,i] = array_list[i]       
    return mat

def adaptive_file_convert(adaptive_file,blen):
    """
    Convert adaptive-K file to P x P matrix where P is the number of SNPs
    Arguments:
    adaptive_file : Npz file containing adaptive array and indices of array endpoints
    blen : Block length used in the function to generate adaptive arrays
    """
    
    # Load in concatenated adaptive array and indices of array endpoints
    inter = np.load(adaptive_file, allow_pickle=True)
    adaptive_mat = inter['adaptive_ld_mat']
    idxs = inter['idx']
    
    array_list = []
    
    # Create list of arrays using endpoints from idxs
    for i in tqdm(range(idxs.size-1)):
        array_list.append(adaptive_ld_mat[idxs[i]:idxs[i+1]])
    
    # Convert list of arrays to P x P matrix
    return ragged_convert(array_list,blen)
 

def stack_ragged(array_list, axis=0):
    """
    Method to stack ragged arrays while retaining
    the indices
    """
    lengths = [np.shape(a)[axis] for a in array_list]
    idx = np.cumsum(lengths[:-1])
    idx = np.insert(idx,0,0)
    stacked = np.concatenate(array_list, axis=axis)
    return(stacked, idx)
