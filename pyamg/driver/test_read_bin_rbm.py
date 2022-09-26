import sys
import numpy as np
from read_bin import read_bin_rbm

try:
    rbmname_in  = sys.argv[1]
except:
    print('USAGE: python3 test_read_bin_rbm.py rbm_in')
    sys.exit()

# Load the binary matrix
RBM = read_bin_rbm(rbmname_in)
print(RBM)
