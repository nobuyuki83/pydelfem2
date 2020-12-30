/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "../py_funcs.h"

#include "delfem2/lsmats.h"
#include "delfem2/lsitrsol.h"
#include "delfem2/lsvecx.h"
#include "delfem2/vecxitrsol.h"
#include "delfem2/lsilu_mats.h"
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// ----------------------------

void MatrixSquareSparse_SetPattern(
    dfm2::CMatrixSparse<double>& mss,
    const py::array_t<unsigned int>& psup_ind,
    const py::array_t<unsigned int>& psup)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( psup_ind.ndim()  == 1 );
  assert( psup.ndim()  == 1 );
  assert( psup_ind.shape()[0] == mss.nrowblk+1 );
  mss.SetPattern(psup_ind.data(), psup_ind.shape()[0],
                 psup.data(),     psup.shape()[0]);
}

void MatrixSquareSparse_SetFixBC(
    dfm2::CMatrixSparse<double>& mss,
    const py::array_t<int>& flagbc)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( flagbc.ndim() == 2 );
  assert( flagbc.shape()[0] == mss.nrowblk );
  assert( flagbc.shape()[1] == mss.nrowdim );
  mss.SetFixedBC(flagbc.data());
}


void PyMatSparse_ScaleBlk_LeftRight(
    dfm2::CMatrixSparse<double>& mss,
    const py::array_t<double>& scale)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( scale.ndim() == 1 );
  assert( scale.shape()[0] == mss.nrowblk );
  MatSparse_ScaleBlk_LeftRight(
      mss,
      scale.data());
}

void PyMatSparse_ScaleBlkLen_LeftRight(
    dfm2::CMatrixSparse<double>& mss,
    const py::array_t<double>& scale)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( scale.ndim() == 2 );
  assert( scale.shape()[0] == mss.nrowblk );
  assert( scale.shape()[1] == mss.nrowdim );
  MatSparse_ScaleBlkLen_LeftRight(mss,
                                  scale.data());
}

void PyMatrixSparse_ScaleBlkLen_LeftRight(
    dfm2::CMatrixSparse<double>& mss,
    const py::array_t<double>& scale,
    bool is_sumndimval)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( scale.ndim() == 2 );
  assert( scale.shape()[0] == mss.nrowblk );
  assert( scale.shape()[1] == mss.nrowdim );
  MatSparse_ScaleBlkLen_LeftRight(mss,
      scale.data());
}

void LinearSystem_SetMasterSlave(
    dfm2::CMatrixSparse<double>& mss,
    py::array_t<double>& np_b,
    const py::array_t<unsigned int>& np_ms)
{
  assert( mss.nrowblk == mss.ncolblk );
  assert( mss.nrowdim == mss.ncoldim );
  assert( dfm2::CheckNumpyArray2D(np_b, mss.nrowblk, mss.nrowdim) );
  assert( dfm2::CheckNumpyArray2D(np_ms, np_b.shape()[0], np_b.shape()[1]) );
  SetMasterSlave(
      mss,
      np_ms.data());
  dfm2::setRHS_MasterSlave(np_b.mutable_data(),
      np_b.shape()[0]*np_b.shape()[1], np_ms.data());
}

std::vector<double> PySolve_PCG(
    py::array_t<double>& vec_b,
    py::array_t<double>& vec_x,
    double conv_ratio, unsigned int iteration,
    const dfm2::CMatrixSparse<double>& mat_A,
    const dfm2::CPreconditionerILU<double>& ilu_A)
{
  //  std::cout << "solve pcg" << std::endl;
  assert( vec_x.size() == vec_b.size() );
  assert( vec_x.size() == mat_A.nrowblk*mat_A.nrowdim );
  const unsigned int N = vec_b.size();
  auto buff_vecb = vec_b.request();
  auto buff_vecx = vec_x.request();
  std::vector<double> tmp0(N);
  std::vector<double> tmp1(N);
  return dfm2::Solve_PCG(
      dfm2::CVecXd((double*)buff_vecb.ptr,N),
      dfm2::CVecXd((double*)buff_vecx.ptr,N),
      dfm2::CVecXd(tmp0.data(),N),
      dfm2::CVecXd(tmp1.data(),N),
      conv_ratio,iteration,
      mat_A,ilu_A);
}

std::vector<double> PySolve_PBiCGStab(
    py::array_t<double>& vec_b,
    py::array_t<double>& vec_x,
    double conv_ratio, unsigned int iteration,
    const dfm2::CMatrixSparse<double>& mat_A,
    const dfm2::CPreconditionerILU<double>& ilu_A)
{
  //  std::cout << "solve pcg" << std::endl;
  auto buff_vecb = vec_b.request();
  auto buff_vecx = vec_x.request();
  return Solve_PBiCGStab((double*)buff_vecb.ptr,
                         (double*)buff_vecx.ptr,
                         conv_ratio,iteration,
                         mat_A,ilu_A);
}


void PyPrecILU_SetPattern_ILUk(
    dfm2::CPreconditionerILU<double>&  mat_ilu,
    const dfm2::CMatrixSparse<double>& mss,
    int nlev_fill)
{
  //  mat_ilu.Initialize_ILU0(mss);
  mat_ilu.Initialize_ILUk(mss,
                          nlev_fill);
}

std::tuple<py::array_t<unsigned int>,py::array_t<unsigned int>>
PyAddMasterSlavePattern(
    const py::array_t<unsigned int>& ms_flag,
    const py::array_t<unsigned int>& np_psup_ind0,
    const py::array_t<unsigned int>& np_psup0)
{
  assert(ms_flag.shape()[0] == np_psup_ind0.shape()[0]-1);
  assert(ms_flag.ndim() == 2 );
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_AddMasterSlavePattern(
      psup_ind, psup,
      ms_flag.data(), ms_flag.shape()[1],
      np_psup_ind0.data(), np_psup_ind0.shape()[0], np_psup0.data());
  py::array_t<unsigned int> np_psup_ind((int)psup_ind.size(),psup_ind.data());
  py::array_t<unsigned int> np_psup((int)psup.size(),psup.data());
  return std::make_tuple(np_psup_ind,np_psup);
}

void PyMasterSlave_DistributeValue(
    py::array_t<double>& val,
    const py::array_t<int>& ms_flag)
{
  double* pVal = (double*)(val.request().ptr);
  const int nDoF = ms_flag.size();
  for(int idof=0;idof<nDoF;++idof){
    int jdof = ms_flag.data()[idof];
    if( jdof == -1 ) continue;
    assert( jdof >= 0 && jdof < nDoF );
    pVal[ idof] = pVal[ jdof];
  }
}


void init_ls(py::module &m){
  py::class_<dfm2::CMatrixSparse<double>>(m,"CppMatrixSparse")
      .def(py::init<>())
      .def("initialize", &dfm2::CMatrixSparse<double>::Initialize)
      .def("set_zero",   &dfm2::CMatrixSparse<double>::setZero)
      .def("add_dia",    &dfm2::CMatrixSparse<double>::AddDia);

  m.def("matrixSquareSparse_setPattern",      &MatrixSquareSparse_SetPattern);
  m.def("matrixSquareSparse_setFixBC",        &MatrixSquareSparse_SetFixBC);
  m.def("cppMatSparse_ScaleBlk_LeftRight",    &PyMatSparse_ScaleBlk_LeftRight);
  m.def("cppMatSparse_ScaleBlkLen_LeftRight", &PyMatSparse_ScaleBlkLen_LeftRight);
  m.def("masterSlave_distributeValue",        &PyMasterSlave_DistributeValue);
  m.def("cppAddMasterSlavePattern",           &PyAddMasterSlavePattern);

  py::class_<dfm2::CPreconditionerILU<double>>(m,"PreconditionerILU")
      .def(py::init<>())
      .def("ilu_decomp", &dfm2::CPreconditionerILU<double>::DoILUDecomp)
      .def("set_value",  &dfm2::CPreconditionerILU<double>::SetValueILU);

  m.def("cppPrecILU_SetPattern_ILUk",    &PyPrecILU_SetPattern_ILUk);

  m.def("linearSystem_setMasterSlave",   &LinearSystem_SetMasterSlave);
  m.def("linsys_solve_pcg",              &PySolve_PCG);
  m.def("linsys_solve_bicgstab",         &PySolve_PBiCGStab);
}
