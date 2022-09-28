"""Support for aggregation-based AMG."""

# Flag to activate Andrea & Carlo DEBUG prints
DEBUG_AC = False

if DEBUG_AC:
    from numpy import savetxt
    from scipy.io import mmwrite

from warnings import warn
import numpy as np
from scipy.sparse import csr_matrix, isspmatrix_csr, isspmatrix_bsr,\
    csc_matrix, SparseEfficiencyWarning, find, bsr_matrix
from scipy.sparse import eye as speye

from ..multilevel import MultilevelSolver
from ..relaxation.smoothing import change_smoothers
from ..relaxation.utils import relaxation_as_linear_operator
from ..util.utils import scale_T, get_Cpt_params, \
    eliminate_diag_dom_nodes, get_blocksize, \
    levelize_strength_or_aggregation, \
    levelize_smooth_or_improve_candidates
from ..strength import classical_strength_of_connection,\
    symmetric_strength_of_connection, evolution_strength_of_connection,\
    energy_based_strength_of_connection, distance_strength_of_connection,\
    algebraic_distance, affinity_distance
from .aggregate import standard_aggregation, naive_aggregation, \
    lloyd_aggregation
from .tentative import fit_candidates
from .smooth import energy_prolongation_smoother
from ..classical import split
from pyamg.amg_core import cptBAMGProl
from .mkPattern import mkPatt
from .EMIN import EMIN

def energymin_cf_solver(A, B=None, BH=None,
                        symmetry='hermitian', strength='symmetric',
                        aggregate='standard', smooth='energy',
                        prolongation='injection',
                        presmoother=('block_gauss_seidel',
                                     {'sweep': 'symmetric'}),
                        postsmoother=('block_gauss_seidel',
                                      {'sweep': 'symmetric'}),
                        improve_candidates=('block_gauss_seidel',
                                            {'sweep': 'symmetric',
                                             'iterations': 4}),
                        max_levels=10, max_coarse=10,
                        diagonal_dominance=False, keep=False, **kwargs):
    """Create a multilevel solver using energy-min AMG

    See the notes below, for the major differences with the classical-style
    smoothed aggregation solver in aggregation.smoothed_aggregation_solver.

    Parameters
    ----------
    A : csr_matrix, bsr_matrix
        Sparse NxN matrix in CSR or BSR format

    B : None, array_like
        Right near-nullspace candidates stored in the columns of an NxK array.
        K must be >= the blocksize of A (see reference [2011OlScTu]_). The default value
        B=None is equivalent to choosing the constant over each block-variable,
        B=np.kron(np.ones((A.shape[0]/get_blocksize(A), 1)), np.eye(get_blocksize(A)))

    BH : None, array_like
        Left near-nullspace candidates stored in the columns of an NxK array.
        BH is only used if symmetry='nonsymmetric'.  K must be >= the
        blocksize of A (see reference [2011OlScTu]_). The default value B=None is
        equivalent to choosing the constant over each block-variable,
        B=np.kron(np.ones((A.shape[0]/get_blocksize(A), 1)), np.eye(get_blocksize(A)))

    symmetry : string
        'symmetric' refers to both real and complex symmetric
        'hermitian' refers to both complex Hermitian and real Hermitian
        'nonsymmetric' i.e. nonsymmetric in a hermitian sense
        Note that for the strictly real case, symmetric and hermitian are
        the same
        Note that this flag does not denote definiteness of the operator.

    strength : list
        Method used to determine the strength of connection between unknowns of
        the linear system.  Method-specific parameters may be passed in using a
        tuple, e.g. strength=('symmetric',{'theta' : 0.25 }). If strength=None,
        all nonzero entries of the matrix are considered strong.

    aggregate : list
        Method used to aggregate nodes.

    smooth : list
        Method used to smooth the tentative prolongator.  Method-specific
        parameters may be passed in using a tuple, e.g.  smooth=
        ('energy',{'krylov' : 'gmres'}).  Only 'energy' and None are valid
        prolongation smoothing options.

    presmoother : tuple, string, list
        Defines the presmoother for the multilevel cycling.  The default block
        Gauss-Seidel option defaults to point-wise Gauss-Seidel, if the matrix
        is CSR or is a BSR matrix with blocksize of 1.  See notes below for
        varying this parameter on a per level basis.

    postsmoother : tuple, string, list
        Same as presmoother, except defines the postsmoother.

    improve_candidates : tuple, string, list
        The ith entry defines the method used to improve the candidates B on
        level i.  If the list is shorter than max_levels, then the last entry
        will define the method for all levels lower.  If tuple or string, then
        this single relaxation descriptor defines improve_candidates on all
        levels.
        The list elements are relaxation descriptors of the form used for
        presmoother and postsmoother.  A value of None implies no action on B.

    max_levels : integer
        Maximum number of levels to be used in the multilevel solver.

    max_coarse : integer
        Maximum number of variables permitted on the coarse grid.

    diagonal_dominance : bool, tuple
        If True (or the first tuple entry is True), then avoid coarsening
        diagonally dominant rows.  The second tuple entry requires a
        dictionary, where the key value 'theta' is used to tune the diagonal
        dominance threshold.

    keep : bool
        Flag to indicate keeping extra operators in the hierarchy for
        diagnostics.  For example, if True, then strength of connection (C),
        tentative prolongation (T), aggregation (AggOp), and arrays
        storing the C-points (Cpts) and F-points (Fpts) are kept at
        each level.

    Other Parameters
    ----------------
    cycle_type : ['V','W','F']
        Structrure of multigrid cycle
    coarse_solver : ['splu', 'lu', 'cholesky, 'pinv', 'gauss_seidel', ... ]
        Solver used at the coarsest level of the MG hierarchy.
        Optionally, may be a tuple (fn, args), where fn is a string such as
        ['splu', 'lu', ...] or a callable function, and args is a dictionary of
        arguments to be passed to fn.

    Returns
    -------
    ml : MultilevelSolver
        Multigrid hierarchy of matrices and prolongation operators

    See Also
    --------
    MultilevelSolver, aggregation.smoothed_aggregation_solver,
    classical.ruge_stuben_solver

    Notes
    -----
         - Here, AMG is "classial" in that it preserves an identity block
           in the interpolation operator, P.  This identity block can represent
           root-nodes of aggregates, or be some other arbitrary CF splitting.

         - Only smooth={'energy', None} is supported for prolongation
           smoothing.

         - The additional parameters are passed through as arguments to
           MultilevelSolver.  Refer to pyamg.MultilevelSolver for additional
           documentation.

         - At each level, four steps are executed in order to define the coarser
           level operator.

           1. Matrix A is given and used to derive a strength matrix, C.

           2. Based on the strength matrix, indices are grouped or aggregated.

           3. The aggregates define coarse nodes and a tentative prolongation
              operator T is defined by injection

           4. The tentative prolongation operator is smoothed by a relaxation
              scheme to improve the quality and extent of interpolation from the
              aggregates to fine nodes.

         - The parameters smooth, strength, aggregate, presmoother, postsmoother
           can be varied on a per level basis.  For different methods on
           different levels, use a list as input so that the i-th entry defines
           the method at the i-th level.  If there are more levels in the
           hierarchy than list entries, the last entry will define the method
           for all levels lower.

           Examples are:
           smooth=[('jacobi', {'omega':1.0}), None, 'jacobi']
           presmoother=[('block_gauss_seidel', {'sweep':symmetric}), 'sor']
           aggregate=['standard', 'naive']
           strength=[('symmetric', {'theta':0.25}), ('symmetric', {'theta':0.08})]

         - Predefined strength of connection and aggregation schemes can be
           specified.  These options are best used together, but aggregation can
           be predefined while strength of connection is not.

           For predefined strength of connection, use a list consisting of
           tuples of the form ('predefined', {'C' : C0}), where C0 is a
           csr_matrix and each degree-of-freedom in C0 represents a supernode.
           For instance to predefine a three-level hierarchy, use
           [('predefined', {'C' : C0}), ('predefined', {'C' : C1}) ].

           Similarly for predefined aggregation, use a list of tuples.  For
           instance to predefine a three-level hierarchy, use [('predefined',
           {'AggOp' : Agg0}), ('predefined', {'AggOp' : Agg1}) ], where the
           dimensions of A, Agg0 and Agg1 are compatible, i.e.  Agg0.shape[1] ==
           A.shape[0] and Agg1.shape[1] == Agg0.shape[0].  Each AggOp is a
           csr_matrix.

           Because this is solver needs a CF splitting, if a member of the predefined
           aggregation list is predefined, it must be of the form
           ('predefined', {'AggOp' : Agg, 'Cnodes' : Cnodes}).

    Examples
    --------
    >>> from pyamg import energymin_cf_solver
    >>> from pyamg.gallery import poisson
    >>> from scipy.sparse.linalg import cg
    >>> import numpy as np
    >>> A = poisson((100, 100), format='csr')           # matrix
    >>> b = np.ones((A.shape[0]))                       # RHS
    >>> ml = energymin_cf_solver(A)                     # AMG solver
    >>> M = ml.aspreconditioner(cycle='V')              # preconditioner
    >>> x, info = cg(A, b, tol=1e-8, maxiter=30, M=M)   # solve with CG

    References
    ----------
    .. [1996VaMa] Vanek, P. and Mandel, J. and Brezina, M.,
       "Algebraic Multigrid by Smoothed Aggregation for
       Second and Fourth Order Elliptic Problems",
       Computing, vol. 56, no. 3, pp. 179--196, 1996.
       http://citeseer.ist.psu.edu/vanek96algebraic.html
    .. [2011OlScTu] Olson, L. and Schroder, J. and Tuminaro, R.,
       "A general interpolation strategy for algebraic
       multigrid using energy minimization", SIAM Journal
       on Scientific Computing (SISC), vol. 33, pp.
       966--991, 2011.

    """
    if not (isspmatrix_csr(A) or isspmatrix_bsr(A)):
        try:
            A = csr_matrix(A)
            warn('Implicit conversion of A to CSR',
                 SparseEfficiencyWarning)
        except BaseException as e:
            raise TypeError('Argument A must have type csr_matrix, '
                            'bsr_matrix, or be convertible to csr_matrix') from e

    A = A.asfptype()

    if symmetry not in ('symmetric', 'hermitian', 'nonsymmetric'):
        raise ValueError('Expected "symmetric", "nonsymmetric" '
                         'or "hermitian" for the symmetry parameter.')
    A.symmetry = symmetry

    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix')
    # Right near nullspace candidates use constant for each variable as default
    if B is None:
        B = np.kron(np.ones((int(A.shape[0]/get_blocksize(A)), 1), dtype=A.dtype),
                    np.eye(get_blocksize(A)))
    else:
        B = np.asarray(B, dtype=A.dtype)
        if len(B.shape) == 1:
            B = B.reshape(-1, 1)
        if B.shape[0] != A.shape[0]:
            raise ValueError('The near null-space modes B have incorrect \
                              dimensions for matrix A')
        if B.shape[1] < get_blocksize(A):
            raise ValueError('B.shape[1] must be >= the blocksize of A')

    # Left near nullspace candidates
    if A.symmetry == 'nonsymmetric':
        if BH is None:
            BH = B.copy()
        else:
            BH = np.asarray(BH, dtype=A.dtype)
            if len(BH.shape) == 1:
                BH = BH.reshape(-1, 1)
            if BH.shape[1] != B.shape[1]:
                raise ValueError('The number of left and right near \
                                  null-space modes B and BH, must be equal')
            if BH.shape[0] != A.shape[0]:
                raise ValueError('The near null-space modes BH have \
                                  incorrect dimensions for matrix A')

    # Levelize the user parameters, so that they become lists describing the
    # desired user option on each level.
    max_levels, max_coarse, strength =\
        levelize_strength_or_aggregation(strength, max_levels, max_coarse)
    max_levels, max_coarse, aggregate =\
        levelize_strength_or_aggregation(aggregate, max_levels, max_coarse)
    improve_candidates =\
        levelize_smooth_or_improve_candidates(improve_candidates, max_levels)
    smooth = levelize_smooth_or_improve_candidates(smooth, max_levels)
    prolongation = levelize_smooth_or_improve_candidates(prolongation, max_levels)

    # Construct multilevel structure
    levels = []
    levels.append(MultilevelSolver.Level())
    levels[-1].A = A          # matrix

    # Append near nullspace candidates
    levels[-1].B = B          # right candidates
    if A.symmetry == 'nonsymmetric':
        levels[-1].BH = BH    # left candidates

    while len(levels) < max_levels and \
            int(levels[-1].A.shape[0]/get_blocksize(levels[-1].A)) > max_coarse:
        _extend_hierarchy(levels, strength, aggregate, smooth,
                          improve_candidates, prolongation, diagonal_dominance, keep)

    ml = MultilevelSolver(levels, **kwargs)
    change_smoothers(ml, presmoother, postsmoother)
    return ml


def _extend_hierarchy(levels, strength, aggregate, smooth, improve_candidates,
                      prolongation, diagonal_dominance=False, keep=True):
    """Extend the multigrid hierarchy.

    Service routine to implement the strength of connection, aggregation,
    tentative prolongation construction, and prolongation smoothing.  Called by
    smoothed_aggregation_solver.

    """
    def unpack_arg(v):
        if isinstance(v, tuple):
            return v[0], v[1]
        return v, {}

    A = levels[-1].A
    B = levels[-1].B
    if A.symmetry == 'nonsymmetric':
        AH = A.H.asformat(A.format)
        BH = levels[-1].BH

    # Compute the strength-of-connection matrix C, where larger
    # C[i, j] denote stronger couplings between i and j.
    fn, kwargs = unpack_arg(strength[len(levels)-1])
    if fn == 'symmetric':
        C = symmetric_strength_of_connection(A, **kwargs)
    elif fn == 'classical':
        C = classical_strength_of_connection(A, **kwargs)
    elif fn == 'distance':
        C = distance_strength_of_connection(A, **kwargs)
    elif fn in ('ode', 'evolution'):
        if 'B' in kwargs:
            C = evolution_strength_of_connection(A, **kwargs)
        else:
            C = evolution_strength_of_connection(A, B, **kwargs)
    elif fn == 'energy_based':
        C = energy_based_strength_of_connection(A, **kwargs)
    elif fn == 'predefined':
        C = kwargs['C'].tocsr()
    elif fn == 'algebraic_distance':
        C = algebraic_distance(A, **kwargs)
    elif fn == 'affinity':
        C = affinity_distance(A, **kwargs)
    elif fn is None:
        C = A.tocsr()
    else:
        raise ValueError(f'Unrecognized strength of connection method: {str(fn)}')

    # Avoid coarsening diagonally dominant rows
    flag, kwargs = unpack_arg(diagonal_dominance)
    if flag:
        C = eliminate_diag_dom_nodes(A, C, **kwargs)

    # Compute the aggregation matrix AggOp (i.e., the nodal coarsening of A).
    # AggOp is a boolean matrix, where the sparsity pattern for the k-th column
    # denotes the fine-grid nodes agglomerated into k-th coarse-grid node.
    fn, kwargs = unpack_arg(aggregate[len(levels)-1])
    if fn == 'standard':
        AggOp, Cnodes = standard_aggregation(C, **kwargs)
    elif fn == 'naive':
        AggOp, Cnodes = naive_aggregation(C, **kwargs)
    elif fn == 'lloyd':
        AggOp, Cnodes = lloyd_aggregation(C, **kwargs)
    elif fn == 'predefined':
        AggOp = kwargs['AggOp'].tocsr()
        Cnodes = kwargs['Cnodes']
    elif fn == 'RS':
        splitting = split.RS(C, **kwargs)
    elif fn == 'PMIS':
        splitting = split.PMIS(C, **kwargs)
    elif fn == 'PMISc':
        splitting = split.PMISc(C, **kwargs)
    elif fn == 'CLJP':
        splitting = split.CLJP(C, **kwargs)
    elif fn == 'CLJPc':
        splitting = split.CLJPc(C, **kwargs)
    elif fn == 'CR':
        # TODO: adaptive C/F splitting based on  abs(T Bc - Bf).  This could be
        # incorporated as a new CR type splitting method, or be implemented around
        # line 1120 in energy_prolongation_smoothing, where filter_operator is used.
        # It may belong better here, at least philosophically.
        splitting = CR(C, **kwargs)
    else:
        raise ValueError(f'Unrecognized aggregation method: {str(fn)}')

    # If using classical splitting, construct needed Rootnode data structures: Cpt_params
    classical_CF = False
    if fn in ['RS', 'PMIS', 'PMISc', 'CLJP', 'CLJPc', 'CR']:
        # TODO: Not sure how to (and if we want) BSR-supernode support for C/F splittings
        if isspmatrix_bsr(A):
            if A.blocksize[0] > 1:
                raise ValueError(f'Currently, BSR A and classical CF splittings are not supported')

        classical_CF = True
        Cpts = (splitting == 1).nonzero()[0]
        Fpts = (splitting == 0).nonzero()[0]
        print('# of nodes:    ',C.shape[0])
        print('# coarse nodes:',Cpts.shape[0])
        I_C = speye(A.shape[0], A.shape[1], format='csr')
        I_F = I_C.copy()
        I_F.data[Cpts] = 0.0
        I_F.eliminate_zeros()
        I_C = I_C - I_F
        I_C.eliminate_zeros()
        # construct P_I as in get_Cpt_params
        indices = Cpts.copy()
        indptr = np.arange(indices.shape[0]+1)
        ncoarse = Cpts.shape[0]
        P_I = csc_matrix((I_C.data.copy(), indices, indptr),
                          shape=(I_C.shape[0], ncoarse))
        # not sure how to handle the BSR case, so just assume it doesn't exist
        P_I = P_I.tobsr((1,1))
        I_C = I_C.tobsr(blocksize=(1, 1))
        I_F = I_F.tobsr(blocksize=(1, 1))
        Cpt_params =  (True, {'P_I': P_I, 'I_F': I_F, 'I_C': I_C, 'Cpts': Cpts, 'Fpts': Fpts})

    # Improve near nullspace candidates by relaxing on A B = 0
    fn, kwargs = unpack_arg(improve_candidates[len(levels)-1])
    if fn is not None:
        b = np.zeros((A.shape[0], 1), dtype=A.dtype)
        B = relaxation_as_linear_operator((fn, kwargs), A, b) * B
        levels[-1].B = B
        if A.symmetry == 'nonsymmetric':
            BH = relaxation_as_linear_operator((fn, kwargs), AH, b) * BH
            levels[-1].BH = BH

    # Compute tentative T
    if classical_CF:
        fn, kwargs = unpack_arg(prolongation[len(levels)-1])
        if fn != 'least_squares':
            # For classical C/F splittings, energy_prolongation_smoother below will
            # enforce T Bc = B
            T = Cpt_params[1]['P_I'].copy()
            if A.symmetry == 'nonsymmetric':
                TH = Cpt_params[1]['P_I'].copy()
        else:
            if A.symmetry == 'nonsymmetric':
                raise ValueError(f'least_squares does not support nonsymmetric matrices')

            # Set default parameters (only for advanced users)
            itmax_vol = 100
            dist_min = 1
            mmax = B.shape[1]
            maxcond = 1.e13
            tol_vol = 0.01
            eps = 1e-10

            # Set default parameters
            verbosity_LS = 0
            dist_max = 6
            maxrownrm = 5.0
            # Extract arguments from kwargs
            if 'verbosity' in kwargs:
                verbosity_LS = kwargs['verbosity']
            if 'dist' in kwargs:
                dist_max = kwargs['dist']
            if 'max_row_norm' in kwargs:
                maxrownrm = kwargs['max_row_norm']

            # Total number of unknowns
            nn_C = C.shape[0]
            iat_C = C.indptr
            ja_C = C.indices
            # Initialize arrays
            nn_I = nn_C
            nc_I = ncoarse
            nt_I = nn_I*mmax
            iat_I = np.empty((nn_I+1,), dtype=np.int32)
            ja_I = np.empty((nt_I,), dtype=np.int32)
            coef_I = np.empty((nt_I,), dtype=float)
            c_mark = np.empty((nn_I,), dtype=np.int32)

            if verbosity_LS > 0:
                print('-------------------------------------------------------')
                print('Computing Adaptive Least Squares tentative prolongation')

            # Create F/C nodes indicator
            fcnodes = np.ones((nn_C,), dtype=np.int32)
            fcnodes[Fpts] = -1

            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if DEBUG_AC:
                if len(levels) == 1:
                    mmwrite('C_X.mtx',C,symmetry='general')
                    savetxt('TV_X.txt',B, header=str(B.shape))
                    savetxt('fc_X.txt',fcnodes, fmt='%3d')
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            # Perform a QR on the candidates (generally it has only little impact)
            from scipy.linalg import qr
            Q, R = qr(B,mode='economic')
            B = Q

            # Compute prolongation
            ierr = cptBAMGProl( len(levels), verbosity_LS, itmax_vol, dist_min, dist_max,
                                mmax, maxcond, maxrownrm, tol_vol, eps, nn_C, iat_C,
                                ja_C, B.shape[1], fcnodes, B.flatten(), nn_I, nt_I,
                                iat_I, ja_I, coef_I, c_mark )
            if (ierr != 0):
                raise ValueError('Error in cptBAMGProl')

            T = csr_matrix((coef_I, ja_I, iat_I), shape=(nn_I, nc_I))
            T = T.tobsr()
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if DEBUG_AC:
                if len(levels) == 1:
                    mmwrite('T_X.mtx',T,symmetry='general')
            #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    else:
        # Compute the tentative prolongator, T, which is a tentative interpolation
        # matrix from the coarse-grid to the fine-grid.  T exactly interpolates
        # B_fine[:, 0:get_blocksize(A)] = T B_coarse[:, 0:get_blocksize(A)].
        T, dummy = fit_candidates(AggOp, B[:, 0:get_blocksize(A)])
        del dummy
        if A.symmetry == 'nonsymmetric':
            TH, dummyH = fit_candidates(AggOp, BH[:, 0:get_blocksize(A)])
            del dummyH

        # Create necessary root node matrices
        Cpt_params = (True, get_Cpt_params(A, Cnodes, AggOp, T))
        T = scale_T(T, Cpt_params[1]['P_I'], Cpt_params[1]['I_F'])
        if A.symmetry == 'nonsymmetric':
            TH = scale_T(TH, Cpt_params[1]['P_I'], Cpt_params[1]['I_F'])

    # Set coarse grid near nullspace modes as injected fine grid near
    # null-space modes
    B = Cpt_params[1]['P_I'].T*levels[-1].B
    if A.symmetry == 'nonsymmetric':
        BH = Cpt_params[1]['P_I'].T*levels[-1].BH

    # Smooth the tentative prolongator, so that it's accuracy is greatly
    # improved for algebraically smooth error.
    fn, kwargs = unpack_arg(smooth[len(levels)-1])
    if fn == 'energy':
        P = energy_prolongation_smoother(A, T, C, B, levels[-1].B,
                                         Cpt_params=Cpt_params,
                                         force_fit_candidates=classical_CF,
                                         **kwargs)
    elif fn == 'EMIN':
        # Set default parameters (only for advanced users)
        condmax_EMIN = 1.e10
        precType = 'jacobi'
        #  Set default parameters
        verbosity_EMIN = 0
        avg_nnzr = 30
        kpow = 1
        itmax_EMIN = 5
        tol_EMIN = 0.01
        # Extract arguments from kwargs
        if 'verbosity' in kwargs:
            verbosity_EMIN = kwargs['verbosity']
        if 'average_nnzr' in kwargs:
            avg_nnzr = kwargs['average_nnzr']
        if 'power_pattern' in kwargs:
            kpow = kwargs['power_pattern']
        if 'itmax' in kwargs:
            itmax_EMIN = kwargs['itmax']
        if 'tol' in kwargs:
            tol_EMIN = kwargs['tol']

        if verbosity_EMIN > 0:
            print('Improving tentative prolongation through energy minimization')

        fcnodes = -np.ones((C.shape[0],), dtype=np.int32)
        fcnodes[Cpts] = range( 0, len(Cpts) )
        pattern = mkPatt(C,T,avg_nnzr,kpow)
        if verbosity_EMIN > 0:
            print('Pattern avg nnzr (including C nodes): {:5.2f}'.format(
                  pattern.nnz/pattern.shape[0]))

        # Remove coaese rows from pattern
        [ii,jj,pp] = find(pattern)
        pos = np.where(fcnodes[ii]<0)[0]
        pattern = csr_matrix((pp[pos], (ii[pos],jj[pos])),shape=pattern.shape)

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if DEBUG_AC:
            mmwrite('C_' + str(len(levels)) + '.mtx',C)
            mmwrite('PATT_' + str(len(levels)) + '.mtx',pattern)
            mmwrite('A_' + str(len(levels)) + '.mtx',A)
            mmwrite('T_' + str(len(levels)) + '.mtx',T)
            savetxt('TV_' + str(len(levels)) + '.txt',levels[-1].B,
                    header=str(levels[-1].B.shape))
            savetxt('fc_' + str(len(levels)) + '.txt',fcnodes, fmt='%3d')
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        P = EMIN(verbosity_EMIN,itmax_EMIN,tol_EMIN,condmax_EMIN,precType,fcnodes,
                 A,T,levels[-1].B,pattern)
        P = P.tobsr()

        if verbosity_EMIN > 0:
            print('------------------------------------------------------------')

        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if DEBUG_AC:
            mmwrite('Pemin_' + str(len(levels)) + '.mtx',P)
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    elif fn is None:
        P = T
    else:
        raise ValueError(f'Unrecognized prolongation smoother method: {str(fn)}')

    # Compute the restriction matrix R, which interpolates from the fine-grid
    # to the coarse-grid.  If A is nonsymmetric, then R must be constructed
    # based on A.H.  Otherwise R = P.H or P.T.
    symmetry = A.symmetry
    if symmetry == 'hermitian':
        R = P.H
    elif symmetry == 'symmetric':
        R = P.T
    elif symmetry == 'nonsymmetric':
        fn, kwargs = unpack_arg(smooth[len(levels)-1])
        if fn == 'energy':
            R = energy_prolongation_smoother(AH, TH, C, BH, levels[-1].BH,
                                             Cpt_params=Cpt_params, **kwargs)
            R = R.H
        elif fn is None:
            R = T.H
        else:
            raise ValueError(f'Unrecognized prolongation smoother method: {str(fn)}')

    if keep:
        levels[-1].C = C                         # strength of connection matrix
        try: levels[-1].AggOp = AggOp            # aggregation operator, not available if using C/F
        except: pass
        levels[-1].T = T                         # tentative prolongator
        levels[-1].Fpts = Cpt_params[1]['Fpts']  # Fpts
        levels[-1].P_I = Cpt_params[1]['P_I']    # Injection operator
        levels[-1].I_F = Cpt_params[1]['I_F']    # Identity on F-pts
        levels[-1].I_C = Cpt_params[1]['I_C']    # Identity on C-pts

    levels[-1].P = P                             # smoothed prolongator
    levels[-1].R = R                             # restriction operator
    levels[-1].Cpts = Cpt_params[1]['Cpts']      # Cpts (i.e., rootnodes)

    levels.append(MultilevelSolver.Level())
    A = R * A * P                                # Galerkin operator
    A.symmetry = symmetry
    levels[-1].A = A
    levels[-1].B = B                             # right near nullspace candidates

    if A.symmetry == 'nonsymmetric':
        levels[-1].BH = BH                       # left near nullspace candidates
