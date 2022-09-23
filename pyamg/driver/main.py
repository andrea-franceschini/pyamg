import numpy
import time
from numpy import ones, ravel, arange, array, dot, abs, zeros_like, zeros, mean, setdiff1d, median

from pyamg import gallery
from pyamg.aggregation import smoothed_aggregation_solver, rootnode_solver, energymin_cf_solver
from pyamg.util.linalg import norm

from scipy.sparse import isspmatrix_csr

from read_bin_CJ import read_bin_csr, read_bin_rbm

##
# Define function to run one test of a solver with initial guess x0 and RHS b

def run_test(solver, b, x0, cycle, krylov, tol, maxiter):
    '''
    Run solver for A x = b, with x0 as initial guess
    Use cycle type 'cycle', krylov acceleration 'krylov' and tolerance 'tol'
    '''

    # Grab matrix from solver structure
    A = solver.levels[0].A
    
    # Print residual norm before solve
    print("Before Solve: ||b-Ax|| = %1.3e"%(norm( ravel(b) - ravel(A*x0) )) )

    def callback(x):
        print( "%15.6e"%( norm( ravel( b )-ravel( A*x ) ) ) )

    # Carry out solve
    residuals = []
    start = time.time()
    x = solver.solve(b, x0=x0, accel=krylov, residuals=residuals, 
                     tol=tol, maxiter=maxiter, cycle=cycle, callback=callback)
    end = time.time(); solve_time = end-start
    cyc_comp = solver.cycle_complexity(cycle=cycle)
    resid_rate = (residuals[-1]/residuals[0])**(1.0/(len(residuals)))
    work_per_doa = (cyc_comp/abs(numpy.log10(resid_rate)))
    asymp_rate = mean( (array(residuals)[1:]/array(residuals)[:-1])[-5:] )

    # Print residual norm after solve
    print("After Solve: ||b-Ax|| = %1.3e"%(norm( ravel(b) - ravel(A*x) )) )
    
    # Print statistics
    print(solver) 
    print("    System size:                      " + str(A.shape) )
    print("    Avg. Resid Reduction:             %1.2f"%resid_rate )
    print("    Asymptotic Resid Rate:            %1.2f"%asymp_rate )
    print("    Iterations:                       %d"%len(residuals) )
    print("    Cycle Complexity:                 %1.2f"%cyc_comp )
    print("    Work per DOA:                     %1.2f"%work_per_doa )
    print("")
    print("    Setup Time:                       %1.2f"%(setup_time) )
    print("    Solve Time:                       %1.2f"%(solve_time) )
 

#-----------------------------------------------------------------------------------------
# Define inputs
matname_in = 'MATRICES/Cubo_35199.csr.npy'
rbmname_in = 'MATRICES/Cubo_35199.RBM.npy'
#matname_in = 'MATRICES/Cubo_591.csr.npy'
#rbmname_in = 'MATRICES/Cubo_591.RBM.npy'

#matname_in = 'MATRICES/Cubo_246389.csr.npy'
#rbmname_in = 'MATRICES/Cubo_246389.RBM.npy'

#matname_in = 'A.npy'
#rbmname_in = 'TV.npy'


# Read matrix and test space from file
A = read_bin_csr(matname_in)
B = read_bin_rbm(rbmname_in)
#print(isspmatrix_csr(A))

#############
#from scipy.io import mmwrite
#mmwrite("A.mtx", A)
#############

#n = 10
#A = gallery.poisson( (n,n) )#, format='csc')
#B = numpy.ones((A.shape[0],))
BH = B.copy()       # left near nullspace
#A = A.tocsr()

numpy.random.seed(10)
#b = numpy.random.rand(A.shape[0],1)
b = numpy.ones((A.shape[0],))
x0 = zeros_like(b)

##
# Solver Parameters

max_levels=15
#max_coarse=150
max_coarse=500
#max_coarse=15
coarse_solver='pinv'
coarse_solver='cholesky'
symmetry='hermitian'

cycle='V'
maxiter=500
krylov= 'cg'
#krylov= 'gmres'

tol=1e-8

##
# Choose Solver Parameters
smooth = ('energy', {'krylov':'cg', 'maxiter':4, 'degree':1, 'weighting':'diagonal' })

strength = ('symmetric', {'theta' : 0.0})
#strength = ('classical', {'theta' : 0.35})
#strength = ('evolution', {'k':2, 'epsilon':2.0, 'symmetrize_measure':True, 'B':None})

#improve_candidates = [ ('block_gauss_seidel', {'sweep':'symmetric', 'iterations':10}) ]
#improve_candidates = [ ('gauss_seidel', {'sweep':'symmetric', 'iterations':1}) ]
improve_candidates = None

#pre_smoother =('jacobi', {'iterations':1})
#post_smoother=('jacobi', {'iterations':1})
pre_smoother =('sfsai_nsy', {'sweep':'symmetric', 'iterations':1})
post_smoother=('sfsai_nsy', {'sweep':'symmetric', 'iterations':1})
#pre_smoother =('gauss_seidel', {'sweep':'forward', 'iterations':1})
#post_smoother=('gauss_seidel', {'sweep':'backward', 'iterations':1})
#pre_smoother =('block_gauss_seidel', {'sweep':'forward', 'iterations':1})
#post_smoother=('block_gauss_seidel', {'sweep':'backward', 'iterations':1})

CF_for_enmin = 'PMIS'
#CF_for_enmin = 'RS'
#CF_for_enmin = 'standard'

CF_for_RN = 'standard'

##
# Run two tests
dom = (True,{'theta':1.1})

print("\n----------   energymin_cf_solver ----------")
start = time.time()
enmin_cf = energymin_cf_solver(A, B=B, BH=BH, strength=strength, smooth=smooth,
             improve_candidates=improve_candidates, aggregate=CF_for_enmin,
             presmoother=pre_smoother, postsmoother=post_smoother, keep=True,
             max_levels=max_levels, max_coarse=max_coarse,
             coarse_solver=coarse_solver, symmetry=symmetry, diagonal_dominance=dom,
             opts={'BAMG','EMIN_AC'})
end = time.time(); setup_time = end-start
        
run_test(enmin_cf, b, x0, cycle, krylov, tol, maxiter)

#print("\n----------   rootnode_solver ----------")
#start = time.time()
#RN = rootnode_solver(A, B=B, BH=BH, strength=strength, smooth=smooth,
#             improve_candidates=improve_candidates, aggregate=CF_for_RN,
#             presmoother=pre_smoother, postsmoother=post_smoother, keep=True,
#             max_levels=max_levels, max_coarse=max_coarse,
#             coarse_solver=coarse_solver, symmetry=symmetry)
#end = time.time(); setup_time = end-start
#        
#run_test(RN, b, x0, cycle, krylov, tol, maxiter)
