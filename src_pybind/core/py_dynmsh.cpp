/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "../py_funcs.h"

#include "delfem2/mshtopoio.h"

#include "delfem2/dtri2_v2dtri.h"
#include "delfem2/dtri3_v3dtri.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// -------------------------------------------------------------

void
PyMeshDynTri3D_Initialize(
    dfm2::CMeshDynTri3D& mesh,
    const py::array_t<double>& po,
    const py::array_t<unsigned int>& tri)
{
  assert( CheckNumpyArray2D(po, -1, -1) ); // 2d or 3d
  assert( CheckNumpyArray2D(tri, -1, 3) );
  mesh.Initialize(po.data(), po.shape()[0], po.shape()[1],
                  tri.data(), tri.shape()[0]);
}

void
PyMeshDynTri2D_Initialize(
    dfm2::CMeshDynTri2D& mesh,
    const py::array_t<double>& po,
    const py::array_t<unsigned int>& tri)
{
  assert( CheckNumpyArray2D(po, -1, 2) );
  assert( CheckNumpyArray2D(tri, -1, 3) );
  mesh.Initialize(po.data(), po.shape()[0],
                  tri.data(), tri.shape()[0]);
}

void
PySetXY_MeshDynTri2D(
    dfm2::CMeshDynTri2D& mesh,
    const py::array_t<double>& npXY)
{
  assert( CheckNumpyArray2D(npXY, -1, 2) );
  assert(npXY.shape()[1]==2);
  mesh.setXY(npXY.data(), npXY.shape()[0]);
}

void
PyCopyMeshDynTri2D(
    py::array_t<double>& npPos,
    py::array_t<unsigned int>& npElm,
    const dfm2::CMeshDynTri2D& mesh)
{
  assert( CheckNumpyArray2D(npPos, -1, 2) );
  assert( CheckNumpyArray2D(npElm, -1, 3) );
  const unsigned int np = mesh.aEPo.size();
  const unsigned int ntri = mesh.aETri.size();
  assert(npPos.shape()[0]==np);
  assert(npElm.shape()[0]==ntri);
  {
    double* pP = (double*)(npPos.request().ptr);
    for(unsigned int ip=0;ip<np;++ip){
      pP[ip*2+0] = mesh.aVec2[ip].x();
      pP[ip*2+1] = mesh.aVec2[ip].y();
    }
  }
  {
    unsigned int* pT = (unsigned int*)(npElm.request().ptr);
    for(unsigned int it=0;it<ntri;++it){
      pT[it*3+0] = mesh.aETri[it].v[0];
      pT[it*3+1] = mesh.aETri[it].v[1];
      pT[it*3+2] = mesh.aETri[it].v[2];
    }
  }
}

// ----------------------------------------------------------------------

/*
std::tuple<std::vector<double>,std::vector<int>>
PyTriangulation
(const std::vector< std::vector<double> >& aaXY,
 double edge_length)
{
  std::vector<int> loopIP_ind,loopIP;
  std::vector<CVector2> aVec2;
  {
    JArray_FromVecVec_XY(loopIP_ind,loopIP, aVec2,
                         aaXY);
    if( !CheckInputBoundaryForTriangulation(loopIP_ind,aVec2) ){
      std::vector<double> aXY;
      std::vector<int> aTri;
      return std::forward_as_tuple(aXY,aTri);
    }
    FixLoopOrientation(loopIP,
                       loopIP_ind,aVec2);
    if( edge_length > 10e-10 ){
      ResamplingLoop(loopIP_ind,loopIP,aVec2,
                     edge_length );
    }
  }
  std::vector<CEPo2> aPo2D;
  std::vector<ETri> aETri;
  Meshing_SingleConnectedShape2D(aPo2D, aVec2, aETri,
                                 loopIP_ind, loopIP);
  if( edge_length > 1.0e-10 ){
    CInputTriangulation_Uniform param(1.0);
    std::vector<int> aFlgPnt(aVec2.size(),0);  // TODO: make this flag
    std::vector<int> aFlgTri(aETri.size(),0);
    MeshingInside(aPo2D,aETri,aVec2, aFlgPnt,aFlgTri,
                  aVec2.size(), 0, edge_length, param);
  }
  std::vector<double> aXY;
  std::vector<int> aTri;
  MeshTri2D_Export(aXY,aTri, aVec2,aETri);
  return std::forward_as_tuple(aXY,aTri);
}
 */

void
PyMapValue(
    py::array_t<double>& npV,
    dfm2::CCmdRefineMesh& mpr)
{
  /*
  assert(npOut.shape()[0]==(int)mpr.iv_ind.size()-1);
  assert(npIn.shape()[0]==(int)mpr.nv_in);
  assert(npIn.shape()[1]==npIn.shape()[1]);
  const int ndimval = npIn.shape()[1];
  const int np1 = npOut.shape()[0];
  double* pV = (double*)(npOut.request().ptr);
  for(int i=0;i<np1*ndimval;++i){ pV[i] = 0.0; }
  for(int ip=0;ip<np1;++ip){
    for(int iip=mpr.iv_ind[ip];iip<mpr.iv_ind[ip+1];++iip){
      int jp = mpr.iv[iip];
      double w0 = mpr.w[iip];
      for(int idim=0;idim<ndimval;++idim){
        pV[ip*ndimval+idim] += w0*npIn.at(jp,idim);
      }
    }
  }
   */
  assert( npV.ndim() == 2 );
  const int np = npV.shape()[0];
  const int ndim = npV.shape()[1];
  double* pV = (double*)(npV.request().ptr);
  mpr.Interpolate(pV, np, ndim);
}

// --------------------------------------

void init_dynmsh(py::module &m){

  py::class_<dfm2::CMeshDynTri3D>(m, "CppMeshDynTri3D")
  .def(py::init<>())
  .def("minmax_xyz",            &dfm2::CMeshDynTri3D::MinMax_XYZ)
  //
  .def("check",                 &dfm2::CMeshDynTri3D::Check)
  .def("ntri",                  &dfm2::CMeshDynTri3D::nTri)
  .def("delete_tri_edge",       &dfm2::CMeshDynTri3D::DeleteTriEdge)
  .def("insert_point_elem",     &dfm2::CMeshDynTri3D::insertPointElem)
  .def("delaunay_around_point", &dfm2::CMeshDynTri3D::DelaunayAroundPoint);
  
  py::class_<dfm2::CMeshDynTri2D>(m, "CppMeshDynTri2D")
  .def(py::init<>())
  .def("minmax_xyz",            &dfm2::CMeshDynTri2D::MinMax_XYZ)
  //
  .def("check",                 &dfm2::CMeshDynTri2D::Check)
  .def("ntri",                  &dfm2::CMeshDynTri2D::nTri)
  .def("npoint",                &dfm2::CMeshDynTri2D::nPoint)
  .def("delete_tri_edge",       &dfm2::CMeshDynTri2D::DeleteTriEdge)
  .def("insert_point_elem",     &dfm2::CMeshDynTri2D::insertPointElem)
  .def("delaunay_around_point", &dfm2::CMeshDynTri2D::DelaunayAroundPoint)
  .def("meshing_loops",         &dfm2::CMeshDynTri2D::meshing_loops)
  .def("refinementPlan_EdgeLongerThan_InsideCircle",   &dfm2::CMeshDynTri2D::RefinementPlan_EdgeLongerThan_InsideCircle);

  py::class_<dfm2::CCmdRefineMesh>(m, "CppMapper")
  .def(py::init<>());

  m.def("map_value",    &PyMapValue);

  // dyntri
  m.def("meshdyntri3d_initialize",&PyMeshDynTri3D_Initialize);
  m.def("meshdyntri2d_initialize",&PyMeshDynTri2D_Initialize);
  m.def("copyMeshDynTri2D",       &PyCopyMeshDynTri2D);
  m.def("setXY_MeshDynTri2D",     &PySetXY_MeshDynTri2D);
}
