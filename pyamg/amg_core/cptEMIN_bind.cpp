// DO NOT EDIT: this file is generated

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include "cptEMIN.h"

namespace py = pybind11;

template <class I, class R>
int _cptEMIN(
              I verbosity,
                  I itmax,
                    R tol,
                R condmax,
               I precType,
                     I nn,
                     I nc,
                    I ntv,
 py::array_t<I> & fcnodes,
   py::array_t<I> & iat_A,
    py::array_t<I> & ja_A,
  py::array_t<R> & coef_A,
  py::array_t<I> & iat_P0,
   py::array_t<I> & ja_P0,
 py::array_t<R> & coef_P0,
 py::array_t<R> & TV_buff,
py::array_t<I> & iat_patt,
 py::array_t<I> & ja_patt,
py::array_t<I> & iat_Pout,
 py::array_t<I> & ja_Pout,
py::array_t<R> & coef_Pout,
    py::array_t<R> & info
             )
{
    auto py_fcnodes = fcnodes.unchecked();
    auto py_iat_A = iat_A.unchecked();
    auto py_ja_A = ja_A.unchecked();
    auto py_coef_A = coef_A.unchecked();
    auto py_iat_P0 = iat_P0.unchecked();
    auto py_ja_P0 = ja_P0.unchecked();
    auto py_coef_P0 = coef_P0.unchecked();
    auto py_TV_buff = TV_buff.unchecked();
    auto py_iat_patt = iat_patt.unchecked();
    auto py_ja_patt = ja_patt.unchecked();
    auto py_iat_Pout = iat_Pout.mutable_unchecked();
    auto py_ja_Pout = ja_Pout.mutable_unchecked();
    auto py_coef_Pout = coef_Pout.mutable_unchecked();
    auto py_info = info.mutable_unchecked();
    const I *_fcnodes = py_fcnodes.data();
    const I *_iat_A = py_iat_A.data();
    const I *_ja_A = py_ja_A.data();
    const R *_coef_A = py_coef_A.data();
    const I *_iat_P0 = py_iat_P0.data();
    const I *_ja_P0 = py_ja_P0.data();
    const R *_coef_P0 = py_coef_P0.data();
    const R *_TV_buff = py_TV_buff.data();
    const I *_iat_patt = py_iat_patt.data();
    const I *_ja_patt = py_ja_patt.data();
    I *_iat_Pout = py_iat_Pout.mutable_data();
    I *_ja_Pout = py_ja_Pout.mutable_data();
    R *_coef_Pout = py_coef_Pout.mutable_data();
    R *_info = py_info.mutable_data();

    return cptEMIN <I, R>(
                verbosity,
                    itmax,
                      tol,
                  condmax,
                 precType,
                       nn,
                       nc,
                      ntv,
                 _fcnodes, fcnodes.shape(0),
                   _iat_A, iat_A.shape(0),
                    _ja_A, ja_A.shape(0),
                  _coef_A, coef_A.shape(0),
                  _iat_P0, iat_P0.shape(0),
                   _ja_P0, ja_P0.shape(0),
                 _coef_P0, coef_P0.shape(0),
                 _TV_buff, TV_buff.shape(0),
                _iat_patt, iat_patt.shape(0),
                 _ja_patt, ja_patt.shape(0),
                _iat_Pout, iat_Pout.shape(0),
                 _ja_Pout, ja_Pout.shape(0),
               _coef_Pout, coef_Pout.shape(0),
                    _info, info.shape(0)
                          );
}

PYBIND11_MODULE(cptEMIN, m) {
    m.doc() = R"pbdoc(
    Pybind11 bindings for cptEMIN.h

    Methods
    -------
    cptEMIN
    )pbdoc";

    py::options options;
    options.disable_function_signatures();

    m.def("cptEMIN", &_cptEMIN<int, double>,
        py::arg("verbosity"), py::arg("itmax"), py::arg("tol"), py::arg("condmax"), py::arg("precType"), py::arg("nn"), py::arg("nc"), py::arg("ntv"), py::arg("fcnodes").noconvert(), py::arg("iat_A").noconvert(), py::arg("ja_A").noconvert(), py::arg("coef_A").noconvert(), py::arg("iat_P0").noconvert(), py::arg("ja_P0").noconvert(), py::arg("coef_P0").noconvert(), py::arg("TV_buff").noconvert(), py::arg("iat_patt").noconvert(), py::arg("ja_patt").noconvert(), py::arg("iat_Pout").noconvert(), py::arg("ja_Pout").noconvert(), py::arg("coef_Pout").noconvert(), py::arg("info").noconvert(),
R"pbdoc(
MAIN PROGRAM)pbdoc");

}

