// DO NOT EDIT: this file is generated

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "cptBAMGProl.h"

namespace py = pybind11;

template <class I, class R>
int _cptBAMGProl(
            I const level,
        I const verbosity,
        I const itmax_vol,
         I const dist_min,
         I const dist_max,
             I const mmax,
          R const maxcond,
        R const maxrownrm,
          R const tol_vol,
              R const eps,
             I const nn_S,
   py::array_t<I> & iat_S,
    py::array_t<I> & ja_S,
              I const ntv,
 py::array_t<I> & fcnodes,
      py::array_t<R> & TV,
             I const nn_I,
                   I nt_I,
   py::array_t<I> & iat_I,
    py::array_t<I> & ja_I,
  py::array_t<R> & coef_I,
  py::array_t<I> & c_mark
                 )
{
    auto py_iat_S = iat_S.unchecked();
    auto py_ja_S = ja_S.unchecked();
    auto py_fcnodes = fcnodes.unchecked();
    auto py_TV = TV.unchecked();
    auto py_iat_I = iat_I.mutable_unchecked();
    auto py_ja_I = ja_I.mutable_unchecked();
    auto py_coef_I = coef_I.mutable_unchecked();
    auto py_c_mark = c_mark.mutable_unchecked();
    const I *_iat_S = py_iat_S.data();
    const I *_ja_S = py_ja_S.data();
    const I *_fcnodes = py_fcnodes.data();
    const R *_TV = py_TV.data();
    I *_iat_I = py_iat_I.mutable_data();
    I *_ja_I = py_ja_I.mutable_data();
    R *_coef_I = py_coef_I.mutable_data();
    I *_c_mark = py_c_mark.mutable_data();

    return cptBAMGProl <I, R>(
                    level,
                verbosity,
                itmax_vol,
                 dist_min,
                 dist_max,
                     mmax,
                  maxcond,
                maxrownrm,
                  tol_vol,
                      eps,
                     nn_S,
                   _iat_S, iat_S.shape(0),
                    _ja_S, ja_S.shape(0),
                      ntv,
                 _fcnodes, fcnodes.shape(0),
                      _TV, TV.shape(0),
                     nn_I,
                     nt_I,
                   _iat_I, iat_I.shape(0),
                    _ja_I, ja_I.shape(0),
                  _coef_I, coef_I.shape(0),
                  _c_mark, c_mark.shape(0)
                              );
}

PYBIND11_MODULE(cptBAMGProl, m) {
    m.doc() = R"pbdoc(
    Pybind11 bindings for cptBAMGProl.h

    Methods
    -------
    cptBAMGProl
    )pbdoc";

    py::options options;
    options.disable_function_signatures();

    m.def("cptBAMGProl", &_cptBAMGProl<int, double>,
        py::arg("level"), py::arg("verbosity"), py::arg("itmax_vol"), py::arg("dist_min"), py::arg("dist_max"), py::arg("mmax"), py::arg("maxcond"), py::arg("maxrownrm"), py::arg("tol_vol"), py::arg("eps"), py::arg("nn_S"), py::arg("iat_S").noconvert(), py::arg("ja_S").noconvert(), py::arg("ntv"), py::arg("fcnodes").noconvert(), py::arg("TV").noconvert(), py::arg("nn_I"), py::arg("nt_I"), py::arg("iat_I").noconvert(), py::arg("ja_I").noconvert(), py::arg("coef_I").noconvert(), py::arg("c_mark").noconvert(),
R"pbdoc(
MAIN PROGRAM)pbdoc");

}

