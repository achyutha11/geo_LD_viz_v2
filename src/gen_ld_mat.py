#!python 

"""
    Functions to compute the K x P LD matrix band and other LD summary
    statistics
"""
import numpy as np
from tqdm import tqdm

def est_kxp_mat(gt_mat, pop_vec, pop, K=500):
    """
        Arguments:
            gt_mat: a 0/1/2 genotype matrix with P rows and N columns
            pop_vec:  a N-length vector of population labels
            pop: a population label for individuals in the broader dataset
            (e.g. 'CEU')
            K: the number of adjacent snps to consider
    """
    P,N = gt_mat.shape
    assert(pop_vec.size==N)
    kxp_ld_mat = np.zeros(shape=(K,P),dtype=np.float32)
    gt_pop_filt=gt_mat[:,(pop_vec==pop)]
    for i in tqdm(range(P)):
        cur_corr = np.ma.corrcoef(gt_pop_filt[i,:],gt_pop_filt[i:(i+K),:])[0,:-1]
        kxp_ld_mat[:cur_corr.shape[0],i] = cur_corr
    # Squaring the matrix to gen r^2
    kxp_ld_mat = kxp_ld_mat**2
    return(kxp_ld_mat)



# TODO : can we write methods here to take in the output of the above function
# and then generate LD-summaries?


#def mat_split(gt_mat, pop_vec, pop,X):
#     """
#        Arguments:
#            gt_mat: a 0/1/2 genotype matrix with P rows and N columns
#            pop_vec:  a N-length vector of population labels
#            pop: a population label for individuals in the broader dataset
#            (e.g. 'CEU')
#            X: number of segments to split the genotype matrix into
#    """
#    P,N = gt_mat.shape
#    assert(pop_vec.size==N)
#    gt_pop_filt=gt_mat[:,(pop_vec==pop)]
#    for i in tqdm(range(X)):
#        return gt_pop_filt((((P/X)*i):((P/X)*(i+1))),:)
     


