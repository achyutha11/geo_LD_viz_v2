#First attempt at a function to return array of LD scores

from tqdm import tqdm

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
