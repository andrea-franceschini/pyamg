import sys
import numpy as np
from scipy.io import mmread

# Utility to read an mtx matrix and dump it as npy file
try:
    matname_in  = sys.argv[1]
    matname_out = sys.argv[2]
except:
    print('USAGE: python3 mtx2binary.py mat_in mat_out')
    sys.exit()

# Load the mtx matrix
A = mmread(matname_in).tocsr()

# Extract arrays
coef = A.data
ja   = A.indices
iat  = A.indptr

# Dump on file
of = open(matname_out,'wb')
np.save(of,coef)
np.save(of,ja)
np.save(of,iat)
of.close()
