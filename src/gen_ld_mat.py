#!python 

"""
    Functions to compute the K x P LD matrix band and other LD summary
    statistics
    NOTE : we can speed this function up 
"""
import numpy as np
from tqdm import tqdm


def est_kxp_mat(gt_mat, pop_vec, pop, K=500):
    """
        Arguments:
            gt_mat: a 0/1/2 genotype matrix with P rows and N columns (missing is np.nan)
            pop_vec:  a N-length vector of population labels
            pop: a population label for individuals in the broader dataset
            (e.g. 'CEU')
            K: the number of adjacent snps to consider
    """
    P,N = gt_mat.shape
    assert(pop_vec.size==N)
    kxp_ld_mat = np.zeros(shape=(K,P),dtype=np.float32)
    gt_pop_filt=gt_mat[:,(pop_vec==pop)]
    
    # Calculating the alternative allele frequency
    alt_ac = np.nansum(gt_pop_filt, axis=1)
    alt_af = alt_ac / (2*gt_pop_filt.shape[1])
    
    for i in tqdm(range(P)):
        cur_corr = np.ma.corrcoef(gt_pop_filt[i,:],gt_pop_filt[i:(i+K),:])[0,:-1]
        kxp_ld_mat[:cur_corr.shape[0],i] = cur_corr
    
    # Squaring the matrix to gen r^2
    kxp_ld_mat = kxp_ld_mat**2
    return(kxp_ld_mat, alt_af)


def running_mean(x, N):
  """
    Calculate a running mean 
  """
  cumsum = np.nancumsum(x) 
  return (cumsum[N:] - cumsum[:-N]) / float(N)

def adaptive_ld_mat_snp(gt, idx, eps=0.05, n=25, blen=100):
  """
    Function to generate an LD-vector for a particular snp
    Arguments:
     gt_mat = genotype matrix (nsnp x nindiv)
     idx : snp index
     eps : minimum of running mean
     n: number of entries over which to calculate the running mean
     blen : block length of entries to calculate running mean 
  """
  assert(n < blen)
  nsnp, nindiv = gt.shape
  assert(idx < nsnp)
  focal_geno_vec = gt[idx,:]
  x = idx + 1
  r2_vec_tot = []
  going = True
  while going and x < nsnp:
    r2_mat = np.corrcoef(focal_geno_vec, gt[x:(x+blen),:])
    r2_vec = r2_mat[0,1:]
    r2_vec2 = r2_vec**2
    cur_mean_r2 = running_mean(r2_vec2, n)
    r2_vec_tot.append(r2_vec2)
    x += blen
    # Check if any entry in the moving average is < eps
    if np.any(cur_mean_r2 <= eps):
      going = False
  if len(r2_vec_tot) > 1:  
    r2_vec_tot = np.hstack(r2_vec_tot)
  else:
    r2_vec_tot = np.array(r2_vec_tot)
  return(r2_vec_tot)
      
  





