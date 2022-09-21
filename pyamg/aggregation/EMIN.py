import numpy as np
from scipy.sparse import csr_matrix
from pyamg.amg_core import cptEMIN

def EMIN(itmax,tol,condmax,precType,fcnodes,A,P0,TV,pattern):
    '''
    Wrapper for the energy minimization
    '''
    # Unpack A
    nn = A.shape[0]
    iat_A  = A.indptr
    ja_A   = A.indices
    coef_A = A.data

    # Unpack P0
    nc = P0.shape[1]
    iat_P0  = P0.indptr
    ja_P0   = P0.indices
    coef_P0 = P0.data

    # Flatten TV
    ntv = TV.shape[1]
    TV_1d = TV.flatten()

    # Unpack pattern (this will be overwritten by the energy min prolongation)
    iat_patt  = pattern.indptr
    ja_patt   = pattern.indices

    # Select preconditioner
    if precType.lower() == 'jacobi':
        precType = 1
    elif precType.lower() == 'seidel':
        precType = 2
    else:
        raise ValueError('Wrong value for precType')

    # Allocate room for the final prolongation
    nt_Pout = pattern.nnz + nc
    iat_Pout = np.empty((nn+1,), dtype=np.int32)
    ja_Pout = np.empty((nt_Pout,), dtype=np.int32)
    coef_Pout = np.empty((nt_Pout,), dtype=np.float64)
    info = np.empty((8,), dtype=np.float64)

    # Call the c function
    ierr = cptEMIN(itmax,tol,condmax,precType,nn,nc,ntv,fcnodes,iat_A,ja_A,coef_A,
                   iat_P0,ja_P0,coef_P0,TV_1d,iat_patt,ja_patt,iat_Pout,ja_Pout,
                   coef_Pout,info)
    if (ierr != 0):
        raise ValueError('Error in cptEMIN')

    # Pack the final prolongation and exit
    Pout = csr_matrix((coef_Pout, ja_Pout, iat_Pout), shape=(nn, nc))

    return Pout
