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

def ragged_convert(array_list):
    """
        Convert the list of arrays to a P x P LD matrix
    """
    p =  len(array_list)
    ld_mat = np.zeros(shape=(p+1,p+1), dtype=np.float32)
    for i in tqdm(range(p)):
        cur_vec = array_list[i]
        ix = i + 1
        l = cur_vec.size
        ld_mat[ix:(ix+l),i] = cur_vec
    return(ld_mat)  
  
def adaptive_file_convert(adaptive_file):
    """
    Convert adaptive-K file to P x P matrix where P is the number of SNPs
    Arguments:
      adaptive_file : Npz file containing adaptive array and indices of array endpoints
    """
    
    # Load in concatenated adaptive array and indices of array endpoints
    inter = np.load(adaptive_file, allow_pickle=True)
    adaptive_mat = inter['adaptive_ld_mat']
    idxs = inter['idx']
    variables = inter['variables']
    idxs = np.insert(idxs,0,0)

    array_list = []
    
    # Create list of arrays using endpoints from idxs
    for i in tqdm(range(idxs.size-1)):
        array_list.append(adaptive_mat[idxs[i]:idxs[i+1]])
    
    # Convert list of arrays to P x P matrix
    return (ragged_convert(array_list),variables)
  
def stack_ragged(array_list, axis=0):
    """
    Method to stack ragged arrays while retaining
    the indices
    """
    lengths = [np.shape(a)[axis] for a in array_list]
    idx = np.cumsum(lengths[:-1])
    stacked = np.concatenate(array_list, axis=axis)
    return(stacked, idx)

def avg_SNP_dist(full_mat_file,eps):
    """
    Method to find the average distance in SNPs for r2 to go below a certain epsilon threshold
    """
    
    loadin = np.load(full_mat_file)
    mat = loadin['ld_mat']
    sample_mat = np.tril(mat)
    
    avg_list = []

    for i in tqdm(range(sample_mat.shape[0])):
        avg_list.append(np.mean(np.diag(sample_mat,-i)))
    
    for i in tqdm(range(len(avg_list))):
        if avg_list[i] < eps:
            return i
            break


