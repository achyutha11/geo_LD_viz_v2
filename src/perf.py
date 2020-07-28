  
#!python 

"""
    Functions to test the performance of K-banded LD matrices and adaptive-K arrays against full LD matrices
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# def fractions(epsilon):
    
#     # Filtering out and counting number of r2 values > eps in each KxP maatrix and full matrix
#     # First row of KxP matrix not counted, all 1s
#     idx = POP_full_mat > epsilon
#     idx100 = POP_K100_kxp[1:,:] > epsilon
#     idx200 = POP_K200_kxp[1:,:] > epsilon
#     idx500 = POP_K500_kxp[1:,:] > epsilon
#     idx1000 = POP_K1000_kxp[1:,:] > epsilon
#     idx2500 = POP_K2500_kxp[1:,:] > epsilon

#     # Full matrix count divided by 2 due to symmetry
#     tot = np.count_nonzero(idx)/2
#     k100 = np.count_nonzero(idx100)
#     k200 = np.count_nonzero(idx200)
#     k500 = np.count_nonzero(idx500)
#     k1000 = np.count_nonzero(idx1000)
#     k2500 = np.count_nonzero(idx2500)

#     # Returns list of fraction of LD above eps covered by each K value 
#     fr = [k100/tot,k200/tot,k500/tot,k1000/tot,k2500/tot]
    
#     return fr

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
  tot_full = np.nansum(full_ld_mat > epsilon)
  tot_kxp = np.nansum(kxp_mat > epsilon)
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
    cur_k, cur_frac = frac_covered_k_banded(full_LD, kxp_mat, epsilon=epsilon)
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
  

  
  
  
# NOTE : we want to avoid the hardcoding of files in there
def perf_test_plot(POP):
    
    # Loading in npz KxP matrices 
    POP_5MB_K100  = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_K100_pos20000000_25000000_pop'+str(POP)+'.npz'
    POP_5MB_K200  = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_K200_pos20000000_25000000_pop'+str(POP)+'.npz'
    POP_5MB_K500  = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_K500_pos20000000_25000000_pop'+str(POP)+'.npz'
    POP_5MB_K1000 = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_K1000_pos20000000_25000000_pop'+str(POP)+'.npz'
    POP_5MB_K2500 = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_K2500_pos20000000_25000000_pop'+str(POP)+'.npz'
    POP_5MB_full  = '/project2/jnovembre/achyutha11/geo_ld_viz_v2/data/performance_tests/chr22_ac5_full_LD_pos20000000_25000000_pop'+str(POP)+'.npz'
    
    POP_K100 = np.load(POP_5MB_K100)
    POP_K200 = np.load(POP_5MB_K200)
    POP_K500 = np.load(POP_5MB_K500)
    POP_K1000 = np.load(POP_5MB_K1000)
    POP_K2500 = np.load(POP_5MB_K2500)
    POP_full = np.load(POP_5MB_full)
    
    global POP_K100_kxp,POP_K200_kxp,POP_K500_kxp,POP_K1000_kxp,POP_K2500_kxp,POP_full_mat
    POP_K100_kxp = POP_K100['ld_mat']
    POP_K200_kxp = POP_K200['ld_mat']
    POP_K500_kxp = POP_K500['ld_mat']
    POP_K1000_kxp = POP_K1000['ld_mat']
    POP_K2500_kxp = POP_K2500['ld_mat']
    POP_full_mat = (POP_full['ld_mat'])

    # Adding various epsilon values 
    e1 = fractions(0.1)
    e2 = fractions(0.25)
    e3 = fractions(0.5)
    e4 = fractions(0.8)
    K = [100,200,500,1000,2500]
    x=[1,1,1,1,1]
  
    # Plotting performance test results and saving figure
    plt.plot(K,e1,label='0.1')
    plt.plot(K,e2,label='0.25')
    plt.plot(K,e3,label='0.5')
    plt.plot(K,e4,label='0.8')
    plt.plot(K,x,linestyle='--')
    plt.legend(title='Epsilon values')
    plt.xlabel('K value')
    plt.ylabel('Percentage of values detected above epsilon')
    plt.title('Performance tests: K-banded LD matrix on chr22, '+str(POP))
    plt.savefig(str(POP)+'_performance_test.png',dpi=300)
    
