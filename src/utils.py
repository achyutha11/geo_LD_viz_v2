#!python 

"""
    Functions to convert K-banded and adaptive-K matrices to P x P matrix (P is number of SNPs)
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def kxp_convert(kxp_file):
  """
  Convert K x P matrix to P x P matrix where P is the number of SNPs
  Argument:
  kxp_file : K x P npz file path
  """
  
  # Load in K x P matrix
  load = np.load(kxp_file)
  kxpmat = load['ld_mat']

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
