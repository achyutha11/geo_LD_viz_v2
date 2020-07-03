#!python 

"""
    Functions to compute the K x P LD matrix band 
"""

import numpy as np
import allel
from tqdm import tqdm
import pandas as pd
import zarr
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
%matplotlib inline

#Read in zarr:
filename = '../../../jnovembre/data/external_public/1kg_phase3/haps/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.zarr'
callset = zarr.open_group(filename, mode='r')

gt = callset['22/calldata/GT']
ac = callset['22/variants/AC']
pos = callset['22/variants/POS']
numalt = callset['22/variants/numalt']

#Convert to numpy arrays
ac_np = ac[:]
numalt_np = numalt[:]
pos_np = pos[:]

true_gt = gt[:,:,0]+gt[:,:,1]

#population vector
popfile = '../../../jnovembre/data/external_public/1kg_phase3/haps/integrated_call_samples_v3.20130502.ALL.panel'
pop_df = pd.read_csv(popfile,sep='\t')
pop_vec = pop_df['pop'].values

#filtering
snp_index_filt = (numalt_np==1) & (ac_np[:,0]>200)
true_gt_ac_filt = true_gt[snp_idx_filt,:]

#Assigning NaN
X = true_gt_ac_filt 
X = X.astype(np.float32)
X[X<0] = np.nan

#Function to compute KxP matrix
def est_kxp_mat(gt_mat,pop_vec,pop,K)
