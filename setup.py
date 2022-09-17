#!/usr/bin/env python
"""PyAMG: Algebraic Multigrid Solvers in Python.

PyAMG is a library of Algebraic Multigrid (AMG)
solvers with a convenient Python interface.

PyAMG features implementations of:

- Ruge-Stuben (RS) or Classical AMG
- AMG based on Smoothed Aggregation (SA)
- Adaptive Smoothed Aggregation (αSA)
- Compatible Relaxation (CR)
- Krylov methods such as CG, GMRES, FGMRES, BiCGStab, MINRES, etc

PyAMG is primarily written in Python with
supporting C++ code for performance critical operations.
"""

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

amg_core_headers = ['evolution_strength',
                    'graph',
                    'krylov',
                    'linalg',
                    'relaxation',
                    'ruge_stuben',
                    'smoothed_aggregation']

ext_modules = [
    Pybind11Extension(f'pyamg.amg_core.{f}',
                      sources=[f'pyamg/amg_core/{f}_bind.cpp'],
                     )
    for f in amg_core_headers]

ext_modules += [
    Pybind11Extension('pyamg.amg_core.sfsai_nsy',
                      sources=['pyamg/amg_core/sfsai_nsy_bind.cpp',
                               'pyamg/amg_core/sfsai_nsy/compress_nsy_sfsai.cpp',
                               'pyamg/amg_core/sfsai_nsy/cpt_nsy_sfsai_coef.cpp',
                               'pyamg/amg_core/sfsai_nsy/fsai_ScalFact.cpp',
                               'pyamg/amg_core/sfsai_nsy/gather_fullsys.cpp',
                               'pyamg/amg_core/sfsai_nsy/iheapsort.cpp',
                               'pyamg/amg_core/sfsai_nsy/merge_row_patt.cpp',
                               'pyamg/amg_core/sfsai_nsy/mk_pattern.cpp',
                               'pyamg/amg_core/sfsai_nsy/power_patt.cpp',
                               'pyamg/amg_core/sfsai_nsy/transp_Patt.cpp'],
                      include_dirs=['pyamg/amg_core/sfsai_nsy'],
                      libraries=['lapacke'],
                     )
    ]

ext_modules += [
    Pybind11Extension('pyamg.amg_core.cptEMIN',
                      sources=['pyamg/amg_core/cptEMIN_bind.cpp',
                               'pyamg/amg_core/cptEMIN/apply_perm.cpp',
                               'pyamg/amg_core/cptEMIN/copy_Prol.cpp',
                               'pyamg/amg_core/cptEMIN/count_rowterms.cpp',
                               'pyamg/amg_core/cptEMIN/cpt_Trace_Acc.cpp',
                               'pyamg/amg_core/cptEMIN/ddot_par.cpp',
                               'pyamg/amg_core/cptEMIN/DEFL_PCG_matfree.cpp',
                               'pyamg/amg_core/cptEMIN/dnrm2_par.cpp',
                               'pyamg/amg_core/cptEMIN/EMIN_ImpProl.cpp',
                               'pyamg/amg_core/cptEMIN/EMIN_matfree.cpp',
                               'pyamg/amg_core/cptEMIN/gather_B_dump.cpp',
                               'pyamg/amg_core/cptEMIN/gather_B_QR.cpp',
                               'pyamg/amg_core/cptEMIN/gather_f.cpp',
                               'pyamg/amg_core/cptEMIN/KP_spmat.cpp',
                               'pyamg/amg_core/cptEMIN/LinvP_spmat.cpp',
                               'pyamg/amg_core/cptEMIN/load_Jacobi.cpp',
                               'pyamg/amg_core/cptEMIN/mkiat_Tglo.cpp',
                               'pyamg/amg_core/cptEMIN/mkiat_Tloc.cpp',
                               'pyamg/amg_core/cptEMIN/mvjcol.cpp',
                               'pyamg/amg_core/cptEMIN/Orth_Q.cpp',
                               'pyamg/amg_core/cptEMIN/print_Q.cpp',
                               'pyamg/amg_core/cptEMIN/Prol_add_Cnodes.cpp',
                               'pyamg/amg_core/cptEMIN/Transp_Patt.cpp',
                               'pyamg/amg_core/cptEMIN/UinvDP_spmat.cpp',
                               'pyamg/amg_core/cptEMIN/wrCSRmat.cpp'],
                      include_dirs=['pyamg/amg_core/cptEMIN'],
                      libraries=['lapacke'],
                      extra_compile_args=['-fopenmp'],
                      extra_link_args=['-fopenmp'],
                     )
    ]

ext_modules += [
    Pybind11Extension('pyamg.amg_core.tests.bind_examples',
                      sources=['pyamg/amg_core/tests/bind_examples_bind.cpp'],
                     )
    ]

setup(
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)
