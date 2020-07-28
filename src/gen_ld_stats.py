#!python 

"""
    Functions to return an array of LD scores for each SNP from a K x P LD matrix 
"""
from tqdm import tqdm
import numpy as np

def ld_score_snp_k(kxp_mat, idx):
  """
    Calculate ld_score based on k-nearest variants for a given snp
  """  
  # Creating empty array to be filled with LD scores, length = number of SNPs    
  localmat = kxp_mat[:,:idx+1]
  # Taking sum of reverse diagonal of matrix
  diag = np.sum(np.fliplr(localmat).diagonal())
  # Calculating sum of ith column values
  col = np.sum(kxp_mat[:,idx])
  # Adding column and diagonal sums
  score = (col + diag - kxp_mat[0,idx]
  return(score)

def ld_score_k(kxp_mat):
  """
    Calculate LD Scores based on the k-nearest variants
    Arguments:
      kxp_mat : a K x P matrix where entries correspond to r^2 metrics
  """
  k,nsnps = kxp_mat.shape
  scores = np.zeros(nsnps, dtype=np.float32)
  for i in tqdm(range(nsnps)):
    scores[i] = ld_score_snp_k(kxp_mat, i)
  return(scores)
           
           
def ld_score_adaptive():
   
           pass
           
           

def ld_score_snp_pos(kxp_mat, pos, idx, win_size=1e6):
  """
    Compute the LD-Score for a given snp index
    Arguments:
      kxp_mat : a K x P matrix where entries correspond to r^2 metrics
      pos : P-length vector 
      idx : integer for snp index
    
  """
  k, nsnps = kxp_mat.shape
  assert(pos.size == nsnps)
  focal_pos = pos[idx]
  # Get the range of positions 
  min_pos, max_pos = (focal_pos - win_size/2.), (focal_pos + win_size/2.)
  idx_win = np.where((pos >= min_pos) & (pos <= max_pos))[0]
  kxp_mat_loc = [:,idx_win]
  pass         
     
                     
  
  
  
