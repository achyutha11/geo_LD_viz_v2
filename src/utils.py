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

  # Define matrix of zeros
  emp = np.zeros((kxpmat.shape[1],kxpmat.shape[1]))
  
  # Iterate through columns of matrix, appending each column from K x P 
  for i in tqdm(range(kxpmat.shape[1])):
      if i < kxpmat.shape[1]-kxpmat.shape[0]:     
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
    # Creating empty matrix
    mat = np.zeros((len(array_list),len(array_list)))
    
    # Iterating list of arrays and appending each list to appropriate position of P x P matrix
    for i in tqdm(range(len(array_list))):
        if i < len(array_list) - blen: 
            if array_list[i].shape[0]>=blen:
                mat[i+1:array_list[i].shape[0]+i+1,i] = array_list[i]
            elif array_list[i].shape[0]==0:
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
