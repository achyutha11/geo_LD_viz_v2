  
#!python 

"""
    Functions to test the performance of K-banded LD matrices and adaptive-K arrays against full LD matrices
"""

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

def fractions(epsilon):
    
    # Filtering out and counting number of r2 values > eps in each KxP maatrix and full matrix
    # First row of KxP matrix not counted, all 1s
    idx = POP_full_mat > epsilon
    idx100 = POP_K100_kxp[1:,:] > epsilon
    idx200 = POP_K200_kxp[1:,:] > epsilon
    idx500 = POP_K500_kxp[1:,:] > epsilon
    idx1000 = POP_K1000_kxp[1:,:] > epsilon
    idx2500 = POP_K2500_kxp[1:,:] > epsilon

    # Full matrix count divided by 2 due to symmetry
    tot = np.count_nonzero(idx)/2
    k100 = np.count_nonzero(idx100)
    k200 = np.count_nonzero(idx200)
    k500 = np.count_nonzero(idx500)
    k1000 = np.count_nonzero(idx1000)
    k2500 = np.count_nonzero(idx2500)

    # Returns list of fraction of LD above eps covered by each K value 
    fr = [k100/tot,k200/tot,k500/tot,k1000/tot,k2500/tot]
    
    return fr
    
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
    
