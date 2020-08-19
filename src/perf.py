  
#!python 

"""
    Functions to test the performance of K-banded LD matrices and adaptive-K arrays against full LD matrices
"""

import numpy as np
from tqdm import tqdm
import copy
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
  
  
    
def frac_cov_k_band_true_r2(R2_true, k=-2, epsilon=1e-2):
  """
    Calculating the fraction of relevant entries kept 
    after setting a band lower K to all 0's
    
    Arguments:
      R2_true: an n x n matrix of the 
      k: the band below which you want to set all entries to be < 0.0
      epsilon: the lower threshold for r2 detection
  """
  n, _ = R2_true.shape
  assert(k < -1)
  assert(n > 3)
  assert((epsilon > 0) & (epsilon <= 1.0))
  x_tril, y_tril = np.tril_indices(n, k=-1)
  x_tril_k, y_tril_k = np.tril_indices(n, k=k)
  R2_x = np.tril(R2_true,k=-1)
  R2_copy = copy.deepcopy(R2_x)
  # Setting the lower elements as 0
  R2_copy[x_tril_k, y_tril_k] = 0.0
  nonzero_r2 = R2_x[R2_x > 0.0]
  nonzero_r2_filt = R2_copy[R2_copy > 0.0]
  frac = np.sum(nonzero_r2_filt > epsilon) / np.sum(nonzero_r2 > epsilon)
  return(frac, nonzero_r2_filt.shape[0])
  

def corrcoef_PxP(R2_true, R2_inf):
  """
    Calculating the correlation coefficient between the lower-triangular entries of an LD matrix and an inferred LD matrix
    Arguments:
      R2_true : n x n sample LD matrix 
      R2_inf : n x n  inferred LD matrix 
  """
  n, _ = R2_true.shape
  n_inf, _ = R2_inf.shape
  # Checking the dimensions and making sure there are enough variants
  assert(n == n_inf)
  assert(n > 3)
  # Get the lower-triangular indices
  x_tril, y_tril = np.tril_indices(n, k=-1)
  r2_true = R2_true[x_tril, y_tril]
  r2_inf = R2_inf[x_tril, y_tril]
  # compute the correlation coefficient to return
  r = np.corrcoef(r2_true, r2_inf)[0,1]
  return(r)

def perf_adaptive(R2_true,converted_adaptive_file,epsilon):
  
  """
    Testing the performance of an adaptive array against the full LD matrix
    Arguments:
    R2_true = n x n sample LD matrix
    converted_adaptive_file = output of adaptive_file_convert on original adaptive file
    epsilon = r2 threshold for performance testing
  """
  assert(len(converted_adaptive_file)==2)
  pop = converted_adaptive_file[1][0]
  region = converted_adaptive_file[1][1]
  eps_adapt = converted_adaptive_file[1][2]
  n_adapt = converted_adaptive_file[1][3]
  blen_adapt = converted_adaptive_file[1][4]
  assert((epsilon > 0) & (epsilon <= 1.0))
  n = R2_true.shape[0]
  mat = converted_adaptive_file[0][:n,:n]
  m = mat.shape[0]
  assert(n==m)
  R2_x = np.tril(R2_true,k=-1)
  nonzero_r2_true = R2_x[R2_x > 0.0]
  nonzero_r2_adaptive = mat[mat > 0.0]
  corrcoef = corrcoef_PxP(R2_true,mat)
  return (pop,region,eps_adapt,n_adapt,blen_adapt,epsilon,np.sum(nonzero_r2_adaptive > epsilon),np.sum(nonzero_r2_true > epsilon),
            np.count_nonzero(nonzero_r2_adaptive),np.count_nonzero(nonzero_r2_true),corrcoef)

def perf_adaptive_v2(R2_true,converted_adaptive_file,epsilon):
  
  """
    Testing the performance of an adaptive array against the full LD matrix (WITHOUT BLEN)
    Arguments:
    R2_true = n x n sample LD matrix
    converted_adaptive_file = output of adaptive_file_convert on original adaptive file
    epsilon = r2 threshold for performance testing
  """
  assert(len(converted_adaptive_file)==2)
  pop = converted_adaptive_file[1][0]
  region = converted_adaptive_file[1][1]
  eps_adapt = converted_adaptive_file[1][2]
  n_adapt = converted_adaptive_file[1][3]
  assert((epsilon > 0) & (epsilon <= 1.0))
  n = R2_true.shape[0]
  mat = converted_adaptive_file[0][:n,:n]
  m = mat.shape[0]
  assert(n==m)
  R2_x = np.tril(R2_true,k=-1)
  nonzero_r2_true = R2_x[R2_x > 0.0]
  nonzero_r2_adaptive = mat[mat > 0.0]
  corrcoef = corrcoef_PxP(R2_true,mat)
  return (pop,region,eps_adapt,n_adapt,epsilon,np.sum(nonzero_r2_adaptive > epsilon),np.sum(nonzero_r2_true > epsilon),
            np.count_nonzero(nonzero_r2_adaptive),np.count_nonzero(nonzero_r2_true),corrcoef)

def n_finder(gt_mat, x, eps):
  """
  Find the number of SNPs from a focal SNP for the mean r2 value to go below a certain epsilon
  Also the mean r2 at each SNP distance from a focal SNP
  """
  numsnps = int(0.95*gt_mat.shape[0])
  assert(x < numsnps)
  assert(eps > 0 and eps<= 1.0)
  indices = np.random.choice(numsnps, size=x, replace=False)
  n = 0
  avg_list = []
  going = True
  while going:
      r2_list = [np.corrcoef(gt_mat[i,:],gt_mat[i+n,:])[0,1]**2 for i in indices]
      avg_list.append(np.mean(r2_list))
      n += 1
      if np.mean(r2_list) < eps:
          going = False 
  return n,avg_list
