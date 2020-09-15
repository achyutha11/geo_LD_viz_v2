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


# -------- Statistics on adaptive-K arrays --------- #

def ld_score_adaptive(npz_file,snp_pos):
    """
    Function to take in particular adaptive file and SNP position and return LD score for that SNP.
    
    Arguments:
    npz_file : Saved compressed file containing adaptive arrays
    snp_pos : Particular position on chromosome, for which LD score is to be calculated    
    """
    
    # Load in adaptive array
    loader = np.load(npz_file)
    adaptive_array = loader['adaptive_ld_mat']
    
    # Split up array into array of arrays, based on idx
    idxs = loader['idx']
    idxs = np.insert(idxs,0,0)
    array_list = []
    
    for i in list(range(idxs.size-1)):
        array_list.append(adaptive_array[idxs[i]:idxs[i+1]])
    
    # Find SNP corresponding to position entered
    positions = loader['positions']
    
    if any(j==snp_pos for i,j in enumerate(positions)):
        for i,j in enumerate(positions):
            if j==snp_pos:
                index = i
    else:
        return 'Invalid position'
        
    # Calculate forward
    forward_ld = sum(array_list[index])
    
    # Calculate backward
    backward_ld = []
    counter = 0
    
    reverse_list = list(range(index))
    
    reverse_list.reverse()
    
    for i in reverse_list:     
        counter+=1      
        if array_list[i].shape[0] <= counter:
            continue
        else:
            backward_ld.append(array_list[i][counter])
        
    backward_sum = sum(backward_ld)
    
    # Total sum, return statement
    total_ld_score = backward_sum + forward_ld
    return total_ld_score
