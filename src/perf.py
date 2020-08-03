  
#!python 

"""
    Functions to test the performance of K-banded LD matrices and adaptive-K arrays against full LD matrices
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def frac_covered_k_banded(full_ld_mat, kxp_mat, epsilon=1e-2):
  """
    Compare the fraction of r2 entries covered 
      between a full LD matrix and a KxP matrix
    Arguments:
      full_ld_mat : p x p numpy.array of r2 values
      kxp : k x p numpy array of r2 values
      epsilon : float value 
    Returns:
      k : value of k 
      frac : fraction of r2 > epsilon detected
  """
  # Check that epsilon is in the appropriate range
  assert((epsilon > 0.) & (epsilon < 1.))
  nsnps,_ = full_ld_mat.shape
  k,p = kxp_mat.shape
  # Check that the number of snps is the same
  assert(nsnps == p)
  # Compute the number of entries 
  tot_full = np.nansum(full_ld_mat > epsilon)/2
  tot_kxp = np.nansum(kxp_mat[2:,:] > epsilon)
  assert(tot_kxp/tot_full <= 1.0)
  
  return(k, tot_kxp/tot_full)

def frac_covered_k_banded_all(full_ld_mat_file, kxp_mat_files, epsilon=1e-2):
  """
    Computing the fraction of r2 values > epsilon covered across many files
    Arguments:
      full_ld_mat_file : a .npz file containing the full LD matrix
      kxp_mat_files : a list of files for the kxp matrices
      epsilon : lower threshold of LD that we want to detect 
    Returns:
      kxs : np.array containing values of k
      fracs : np.array containing values of the fractions
  """
  full_LD = np.load(full_ld_mat_file)
  full_mat = full_LD['ld_mat']
  n = len(kxp_mat_files)
  ks = np.zeros(n)
  fracs = np.zeros(n)
  i = 0
  for file in tqdm(kxp_mat_files):
    kxp_mat_data = np.load(file)
    kxp_mat = kxp_mat_data['ld_mat']
    cur_k, cur_frac = frac_covered_k_banded(full_mat, kxp_mat, epsilon=epsilon)
    ks[i] = cur_k
    fracs[i] = cur_frac
    i += 1
  return(ks, fracs)


def frac_covered_adaptive(full_ld_mat, adaptive_ld_mat, idx_adaptive, epsilon=1e-2):
  """
    Compare the fraction of r2 entries covered 
      between a full LD matrix and an adaptive LD matrix
    Arguments:
      full_ld_mat : p x p numpy.array of r2 values
      kxp : k x p numpy array of r2 values
      epsilon : float value 
  """
  # Check that epsilon is in the appropriate range
  assert((epsilon > 0.) & (epsilon < 1.))
  nsnps,_ = full_ld_mat.shape
  nsnps_inf = idx_adaptive.size
  assert(nsnps == nsnps_inf)
  # Make sure 
  assert(adaptive_ld_mat.size <= (nsnps**2))
  tot_full = np.nansum(full_ld_mat > epsilon)
  tot_adaptive = np.nansum(adaptive_ld_mat > epsilon)
  
  assert(tot_adaptive/tot_full <= 1.0)
  return(tot_adaptive / tot_full)
  

  
 
