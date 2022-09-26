// DO NOT EDIT: this file is generated

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "sfsai_nsy.h"

namespace py = pybind11;

template <class I, class R>
int _sfsai_nsy(
              I verbosity,
                   I kpow,
               I nnzr_max,
               R tau_pref,
               R tau_post,
                   I nn_A,
                   I nt_A,
   py::array_t<I> & iat_A,
    py::array_t<I> & ja_A,
  py::array_t<R> & coef_A,
  py::array_t<I> & iat_FL,
   py::array_t<I> & ja_FL,
 py::array_t<R> & coef_FL,
  py::array_t<I> & iat_FU,
   py::array_t<I> & ja_FU,
 py::array_t<R> & coef_FU
               )
{
    auto py_iat_A = iat_A.unchecked();
    auto py_ja_A = ja_A.unchecked();
    auto py_coef_A = coef_A.unchecked();
    auto py_iat_FL = iat_FL.mutable_unchecked();
    auto py_ja_FL = ja_FL.mutable_unchecked();
    auto py_coef_FL = coef_FL.mutable_unchecked();
    auto py_iat_FU = iat_FU.mutable_unchecked();
    auto py_ja_FU = ja_FU.mutable_unchecked();
    auto py_coef_FU = coef_FU.mutable_unchecked();
    const I *_iat_A = py_iat_A.data();
    const I *_ja_A = py_ja_A.data();
    const R *_coef_A = py_coef_A.data();
    I *_iat_FL = py_iat_FL.mutable_data();
    I *_ja_FL = py_ja_FL.mutable_data();
    R *_coef_FL = py_coef_FL.mutable_data();
    I *_iat_FU = py_iat_FU.mutable_data();
    I *_ja_FU = py_ja_FU.mutable_data();
    R *_coef_FU = py_coef_FU.mutable_data();

    return sfsai_nsy <I, R>(
                verbosity,
                     kpow,
                 nnzr_max,
                 tau_pref,
                 tau_post,
                     nn_A,
                     nt_A,
                   _iat_A, iat_A.shape(0),
                    _ja_A, ja_A.shape(0),
                  _coef_A, coef_A.shape(0),
                  _iat_FL, iat_FL.shape(0),
                   _ja_FL, ja_FL.shape(0),
                 _coef_FL, coef_FL.shape(0),
                  _iat_FU, iat_FU.shape(0),
                   _ja_FU, ja_FU.shape(0),
                 _coef_FU, coef_FU.shape(0)
                            );
}

PYBIND11_MODULE(sfsai_nsy, m) {
    m.doc() = R"pbdoc(
    Pybind11 bindings for sfsai_nsy.h

    Methods
    -------
    sfsai_nsy
    )pbdoc";

    py::options options;
    options.disable_function_signatures();

    m.def("sfsai_nsy", &_sfsai_nsy<int, double>,
        py::arg("verbosity"), py::arg("kpow"), py::arg("nnzr_max"), py::arg("tau_pref"), py::arg("tau_post"), py::arg("nn_A"), py::arg("nt_A"), py::arg("iat_A").noconvert(), py::arg("ja_A").noconvert(), py::arg("coef_A").noconvert(), py::arg("iat_FL").noconvert(), py::arg("ja_FL").noconvert(), py::arg("coef_FL").noconvert(), py::arg("iat_FU").noconvert(), py::arg("ja_FU").noconvert(), py::arg("coef_FU").noconvert(),
R"pbdoc(
MAIN PROGRAM)pbdoc");

}

