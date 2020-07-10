#!python 

"""
    Functions to return an array of LD scores for each SNP from a K x P LD matrix 
"""

from tqdm import tqdm
import numpy as np

# First attempt to create LD score function

def ld_score(LDmat):
    # Creating empty array to be filled with LD scores, length = number of SNPs
    scores = np.zeros(LDmat.shape[1])
    # While loop creates diagonal component of LD score, ends either when first column or last row of LD matrix is reached
    for i in tqdm(range(LDmat.shape[1])):
        col = i 
        row = 0
        diag = 0
        # Adding each cell in the diagonal to diag
        while col>=1 and row<(LDmat.shape[0]-2):
            diag += LDmat[1+row,col-1] 
            col -= 1
            row += 1
        # Adds sum of column to sum of diagonal
        scores[i]=np.sum(LDmat[:,i]) + diag
    return scores

# Second attempt to create LD score function

def ld_score_v2(LDmat):
    # Creating empty array to be filled with LD scores, length = number of SNPs
    scores = np.zeros(LDmat.shape[1])
    # Iterating through each column of K x P matrix
    for i in tqdm(range(LDmat.shape[1])):        
        localmat = LDmat[:,:i+1]
        # Taking sum of reverse diagonal of matrix
        flpdiag = np.fliplr(localmat).diagonal()
        diag = np.sum(flpdiag)
        # Calculating sum of ith column values
        col = np.sum(LDmat[:,i])
        # Adding column and diagonal sums
        scores[i] = col + diag - LDmat[0,i]
    return scores
