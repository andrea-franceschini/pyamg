import sys
import numpy
from numpy import loadtxt

# Utility to read a test space in ascii and dump it in a npy file
try:
    TS_name_in  = sys.argv[1]
    TS_name_out   = sys.argv[2]
except:
    print('USAGE: python TStxt2npy.py TS_name_in TS_name_out')
    sys.exit()

# Read RBM in input
TS = loadtxt(TS_name_in,dtype=numpy.float64)
[ndof,ntv] = TS.shape

# Dump on file
of = open(TS_name_out,'wb')
numpy.save(of,TS)
of.close()
