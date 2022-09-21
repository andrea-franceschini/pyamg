"""Aggregation-based AMG."""

from .adaptive import adaptive_sa_solver
from .aggregate import (standard_aggregation, naive_aggregation,
                        lloyd_aggregation, balanced_lloyd_aggregation)
from .aggregation import smoothed_aggregation_solver
from .tentative import fit_candidates
from .smooth import (jacobi_prolongation_smoother, richardson_prolongation_smoother,
                     energy_prolongation_smoother)
from .rootnode import rootnode_solver
from .energymin_cf_solver import energymin_cf_solver

from .mkPattern import mkPatt
from .EMIN import EMIN

__all__ = ['adaptive_sa_solver',
           'standard_aggregation', 'naive_aggregation',
           'lloyd_aggregation', 'balanced_lloyd_aggregation',
           'smoothed_aggregation_solver',
           'fit_candidates',
           'jacobi_prolongation_smoother', 'richardson_prolongation_smoother',
           'energy_prolongation_smoother',
           'rootnode_solver', 
           'energymin_cf_solver']
