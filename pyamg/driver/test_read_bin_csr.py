import sys
import numpy as np
from read_bin import read_bin_csr

try:
    matname_in  = sys.argv[1]
except:
    print('USAGE: python3 test_read_bin_csr.py mat_in')
    sys.exit()

# Load the binary matrix
A = read_bin_csr(matname_in)
print(A)
