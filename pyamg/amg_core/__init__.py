"""amg_core - a C++ implementation of AMG-related routines."""

from .evolution_strength import (apply_absolute_distance_filter, apply_distance_filter,
                                 min_blocks, evolution_strength_helper,
                                 incomplete_mat_mult_csr)
from .graph import (maximal_independent_set_serial, maximal_independent_set_parallel,
                    vertex_coloring_mis, vertex_coloring_jones_plassmann,
                    vertex_coloring_LDF,
                    cluster_node_incidence, cluster_center,
                    bellman_ford, lloyd_cluster, lloyd_cluster_exact,
                    maximal_independent_set_k_parallel,
                    breadth_first_search, connected_components)

from .krylov import (apply_householders, householder_hornerscheme, apply_givens)
from .linalg import (pinv_array, csc_scale_columns, csc_scale_rows)
from .relaxation import (gauss_seidel, bsr_gauss_seidel, gauss_seidel_indexed,
                         jacobi, bsr_jacobi,
                         jacobi_ne, gauss_seidel_ne, gauss_seidel_nr,
                         block_jacobi, block_gauss_seidel,
                         extract_subblocks, overlapping_schwarz_csr)
from .ruge_stuben import (classical_strength_of_connection_abs,
                          classical_strength_of_connection_min,
                          maximum_row_value,
                          rs_cf_splitting, rs_cf_splitting_pass2,
                          cljp_naive_splitting,
                          rs_direct_interpolation_pass1, rs_direct_interpolation_pass2,
                          cr_helper)
from .smoothed_aggregation import (symmetric_strength_of_connection, standard_aggregation,
                                   naive_aggregation,
                                   fit_candidates,
                                   satisfy_constraints_helper, calc_BtB,
                                   incomplete_mat_mult_bsr, truncate_rows_csr)
from .sfsai_nsy import sfsai_nsy

from .cptEMIN import cptEMIN

__all__ = [
    'apply_absolute_distance_filter',
    'apply_distance_filter',
    'min_blocks',
    'evolution_strength_helper',
    'incomplete_mat_mult_csr',
    #
    'maximal_independent_set_serial',
    'maximal_independent_set_parallel',
    'vertex_coloring_mis',
    'vertex_coloring_jones_plassmann',
    'vertex_coloring_LDF',
    'cluster_node_incidence',
    'cluster_center',
    'bellman_ford',
    'lloyd_cluster',
    'lloyd_cluster_exact',
    'maximal_independent_set_k_parallel',
    'breadth_first_search',
    'connected_components',
    #
    'apply_householders',
    'householder_hornerscheme',
    'apply_givens',
    #
    'pinv_array',
    'csc_scale_columns',
    'csc_scale_rows',
    #
    'gauss_seidel',
    'bsr_gauss_seidel',
    'jacobi',
    'bsr_jacobi',
    'gauss_seidel_indexed',
    'jacobi_ne',
    'gauss_seidel_ne',
    'gauss_seidel_nr',
    'block_jacobi',
    'block_gauss_seidel',
    'extract_subblocks',
    'overlapping_schwarz_csr',
    #
    'classical_strength_of_connection_abs',
    'classical_strength_of_connection_min',
    'maximum_row_value',
    'rs_cf_splitting',
    'rs_cf_splitting_pass2',
    'cljp_naive_splitting',
    'rs_direct_interpolation_pass1',
    'rs_direct_interpolation_pass2',
    'cr_helper',
    #
    'symmetric_strength_of_connection',
    'standard_aggregation',
    'naive_aggregation',
    'fit_candidates',
    'satisfy_constraints_helper',
    'calc_BtB',
    'incomplete_mat_mult_bsr',
    'truncate_rows_csr'
    #
    'sfsai_nsy'
]
