#!python 

"""
    Functions to return an array of LD scores for each SNP from a K x P LD matrix 
"""
import numpy as np
from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr


# -------- Statistics on real genotype matrices --------- #
def ld_score_snp_pos(geno_mat, pos, idx, win_size=1e6):
  """
    Compute the LD-Score for a given snp index
    Args:
      geno_mat : an N x P matrix where entries are the genotypes
      pos : P-length vector 
      idx : integer for snp index
  """
  n, nsnps = geno_mat.shape
  assert(pos.size == nsnps)
  focal_pos = pos[idx]
  # Get the range of positions 
  min_pos, max_pos = (focal_pos - win_size), (focal_pos + win_size)
  idx_win = np.where((pos >= min_pos) & (pos <= max_pos))[0]
  r2_score = np.zeros(idx_win.size)
  i = 0
  for x in idx_win:
    r_test = pearsonr(geno_mat[:,idx], geno_mat[:,x])[0]
    r2_score[i] = (r_test**2)
    i += 1
  r2_score = 1.0 + np.nansum(r2_score)
  return(r2_score)         
        
def ld_score_all_pos(geno_mat, pos, win_size=1e6):
  """
    Compute the LD-Score for all SNPs within a given window
    Args:
      geno_mat : an N x P matrix where entries are the 
      pos: 
  """
  n, nsnps = geno_mat.shape
  assert(pos.size == nsnps)
  ld_scores = np.zeros(nsnps)
  for i in tqdm(range(nsnps)):
    ld_scores[i] = ld_score_snp_pos(geno_mat, pos, idx=i, win_size=win_size)
  return(ld_scores)