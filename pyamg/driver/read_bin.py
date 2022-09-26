import numpy as np
from scipy.sparse import csr_matrix

def read_bin_csr(matname):
    '''
    Function that reads a csr matrix from a binary file
    '''

    # Read numpy arrays from the input file
    of = open(matname,'rb')
    coef = np.load(of)
    ja   = np.load(of)
    iat  = np.load(of)
    of.close()

    # Retrieve dimensions
    nrows = iat.shape[0] - 1
    ncols = max(ja) + 1

    # Create the csr matrix
    A = csr_matrix((coef, ja, iat), shape=(nrows,ncols))

    return A

def read_bin_rbm(rbmname):
    '''
    Function that reads rigid body modes from a binary file
    '''

    # Read numpy arrays from the input file
    of = open(rbmname,'rb')
    RBM = np.load(of)
    of.close()

    return RBM
