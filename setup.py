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
    Pybind11Extension('pyamg.amg_core.tests.bind_examples',
                      sources=['pyamg/amg_core/tests/bind_examples_bind.cpp'],
                     )
    ]

setup(
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)
