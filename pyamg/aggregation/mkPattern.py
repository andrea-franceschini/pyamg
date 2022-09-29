from scipy.sparse import find, csr_matrix
from numpy import argsort

#-----------------------------------------------------------------------------------------

def FilterMat(A,avg_nnzr):
    '''
    Filters the input matrix A leaving an average of avg_nnzr non-zeroes per rows
    '''

    # Extract indices and values
    [ii,jj,aa] = find(A)

    # Select those with maximum absolute value
    p = abs(aa).argsort()
    n = A.shape[0]
    n_retain = n*avg_nnzr
    p = p[-n_retain:]

    # Create a new matrix
    ii = ii[p]
    jj = jj[p]
    aa = aa[p]
    B = csr_matrix((aa,(ii,jj)),shape=A.shape)

    return B

#-----------------------------------------------------------------------------------------

def mkPatt(S,P0_in,avg_nnzr,kpow):
    '''
    Create prolongation pattern using SoC and initial prolongation
    '''

    # Assume S include unitary diagonal

    # Set all ones in P0
    P0 = P0_in.copy()
    P0.data[0:] = 1
    P_out = P0
    nnzr = P_out.nnz / P_out.shape[0]
    print('Power {:2d} nnzr {:10.3f}'.format(0,nnzr))
    for i in range(kpow):
        # Perform the product times S
        P_out = S*P_out
        # Filter P_out
        P_out = FilterMat(P_out,avg_nnzr)
        nnzr = P_out.nnz / P_out.shape[0]
        print('Power {:2d} nnzr {:10.3f}'.format(i+1,nnzr))

    # Transform P_out in a zero pattern
    P_out.data[0:] = 1

    return P_out
