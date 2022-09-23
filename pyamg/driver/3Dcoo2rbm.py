import sys
import numpy
from numpy import loadtxt
from pyamg.util.utils import coord_to_rbm

# Utility to read a set of 3D coordinates and dump corresponding RBM in a npy file
try:
    coord_name_in  = sys.argv[1]
    RBM_name_out   = sys.argv[2]
except:
    print('USAGE: python3 3Dcoo2rbm.py coord_name_in RBM_name_out')
    sys.exit()

# Read coordinates
coord = loadtxt(coord_name_in,dtype=numpy.float64)
[nnodes,ndof] = coord.shape

# Check dimension
if ndof != 3:
    print('This function only supports 3D coordinates')
    sys.exit()

# Create Rigid Body Modes
RBM = coord_to_rbm(nnodes,ndof,coord[:,0],coord[:,1],coord[:,2])

# Dump on file
of = open(RBM_name_out,'wb')
numpy.save(of,RBM)
of.close()
