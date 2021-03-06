/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#ifndef py_funcs_h
#define py_funcs_h

#include <cassert>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace delfem2 {

template<typename T>
bool CheckNumpyArray2D(
    const pybind11::array_t<T> &np0,
    int n0, int n1) {
  if (np0.ndim() != 2) { goto FAIL; }
  if (n0 >= 0) { if (np0.shape()[0] != n0) { goto FAIL; }}
  if (n1 >= 0) { if (np0.shape()[1] != n1) { goto FAIL; }}
  if (np0.strides(0) != (int) sizeof(T) * np0.shape()[1]) { goto FAIL; }
  if (np0.strides(1) != (int) sizeof(T)) { goto FAIL; }
  return true;
  FAIL:
  std::cout << "fail to assert" << std::endl;
  std::cout << "   " << np0.ndim() << std::endl;
  std::cout << "   " << np0.shape()[0] << std::endl;
  std::cout << "   " << np0.shape()[1] << std::endl;
  return false;
}

template<typename T>
bool CheckJArray(
    const pybind11::array_t<T> &psup_ind,
    const pybind11::array_t<T> &psup) {
  const unsigned int np = psup_ind.shape()[0] - 1;
  if (psup_ind.ndim() != 1) { goto FAIL; }
  if (psup.ndim() != 1) { goto FAIL; }
  if (psup_ind.at(0) != 0) { goto FAIL; }
  if (psup.shape()[0] != psup_ind.at(np)) { goto FAIL; }
  if (psup_ind.strides(0) != (int) sizeof(T)) { goto FAIL; }
  if (psup.strides(0) != (int) sizeof(T)) { goto FAIL; }
  return true;
  FAIL:
  return false;
}

}

#endif /* py_funcs_h */
