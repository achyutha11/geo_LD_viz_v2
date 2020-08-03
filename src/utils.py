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
            emp[i+1:i+1+kxpmat[:,i].shape[0],i] = kxpmat[:,i]
        # Accounting for edge cases
        else:
            emp[i+1:,i] = kxpmat[:,i][:kxpmat.shape[1]-i-1]
    
    return emp

def ragged_convert(array_list,blen): 
    """
    Convert list of adaptive-K arrays to P x P matrix where P is the number of SNPs
    Arguments:
    array_list : List of arrays of different lengths
    blen : block length used to generate the list of arrays 
    """
    assert(blen < len(array_list))
    assert(len(array_list) > 1)
    P = len(array_list)
    
    # Creating empty matrix
    mat = np.zeros((P,P))
    
    # Iterating list of arrays and appending each list to appropriate position of P x P matrix
    for i in tqdm(range(P)):
        if i < P - blen: 
            if array_list[i].shape[0] >= blen:
                mat[i+1:array_list[i].shape[0]+i+1,i] = array_list[i]
            elif array_list[i].shape[0] == 0:
                continue
            elif array_list[i].shape[0] > 0 and array_list[i].shape[0] < blen:
                sub_indices = []
                for j in range(array_list[i].shape[0]):
                    sub_indices.append(j)
                if all(array_list[i][m].shape[0]%blen==0 for m in sub_indices):
                    index = array_list[i].shape[0]*blen
                    mat[i+1:index+i+1,i] = np.concatenate(array_list[i]) 
                else:
                    sums = []
                    for k in range(array_list[i].shape[0]):
                        sums.append(array_list[i][k].shape[0])
                    index = sum(sums)
                    mat[i+1:index+i+1,i] = np.concatenate(array_list[i]) 
        else:
            mat[i+1:array_list[i].shape[0]+i+1,i] = array_list[i]       
    return mat

