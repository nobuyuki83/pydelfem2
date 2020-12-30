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
#include "delfem2/jagarray.h"


namespace py = pybind11;
namespace dfm2 = delfem2;

// ----------------------------------------------------------

void PySetTopology_ExtrudeTri2Tet(
    py::array_t<unsigned int>& npTet,
    int nlayer,
    int nXY,
    const py::array_t<unsigned int>& npTri)
{
  assert( dfm2::CheckNumpyArray2D(npTet, -1, 4) );
  assert( dfm2::CheckNumpyArray2D(npTri, -1, 3) );
  dfm2::SetTopology_ExtrudeTri2Tet(npTet.mutable_data(),
      nXY,
      npTri.data(), npTri.shape()[0],
      nlayer);
}

// --------------------------------------------------------

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshTri3D_ReadPly(
    const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Ply(fname, aXYZ, aTri);
  py::array_t<double> np_XYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<unsigned int> np_Tri({(int)aTri.size()/3,3}, aTri.data());
  return std::make_tuple(np_XYZ,np_Tri);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshTri3D_ReadObj(
    const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_Obj3(fname, aXYZ, aTri);
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::make_tuple(npXYZ,npTri);
}

void
PyMeshTri3D_WriteObj(
    const std::string& fname,
    const py::array_t<double>& aXYZ,
    const py::array_t<unsigned int>& aTri)
{
  assert( dfm2::CheckNumpyArray2D(aXYZ, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(aTri, -1, 3) );
  delfem2::Write_Obj(fname,
            aXYZ.data(), aXYZ.shape()[0],
            aTri.data(), aTri.shape()[0]);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshTri3D_ReadNastran(
    const std::string& fname)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  delfem2::Read_MeshTri3D_Nas(aXYZ, aTri, fname.c_str());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  return std::make_tuple(npXYZ,npTri);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshQuad3D_Subviv(
    const py::array_t<double>& aXYZ0,
    const py::array_t<unsigned int>& aQuad0)
{
  assert( dfm2::CheckNumpyArray2D(aXYZ0, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(aQuad0, -1, 4) );
  std::vector<unsigned int> aQuad1;
  std::vector<unsigned int> psupIndQuad0, psupQuad0;
  std::vector<int> aEdgeFace0;
  dfm2::SubdivTopo_MeshQuad(aQuad1,
                            psupIndQuad0,psupQuad0, aEdgeFace0,
                            aQuad0.data(), aQuad0.shape()[0], aXYZ0.shape()[0]);
  //
  std::vector<double> aXYZ1;
  dfm2::SubdivisionPoints_QuadCatmullClark(aXYZ1,
      aQuad1,aEdgeFace0,psupIndQuad0,psupQuad0,
      aQuad0.data(),aQuad0.shape()[0],
      aXYZ0.data(), aXYZ0.shape()[0]);
  py::array_t<double> npXYZ1({(int)aXYZ1.size()/3,3}, aXYZ1.data());
  py::array_t<unsigned int> npQuad1({(int)aQuad1.size()/4,4}, aQuad1.data());
  return std::make_tuple(npXYZ1,npQuad1);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshHex3D_Subviv(
    const py::array_t<double>& aXYZ0,
    const py::array_t<unsigned int>& aHex0)
{
  assert( dfm2::CheckNumpyArray2D(aXYZ0, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(aHex0, -1, 8) );
  std::vector<unsigned int> aHex1;
  std::vector<unsigned int> psupIndHex0, psupHex0;
  std::vector<unsigned int> aQuadHex0;
  dfm2::SubdivTopo_MeshHex(aHex1,
      psupIndHex0, psupHex0,
      aQuadHex0,
      //
      aHex0.data(), aHex0.shape()[0], aXYZ0.shape()[0]);
  std::vector<double> aXYZ1;
  dfm2::SubdivisionPoints_Hex(aXYZ1,
      psupIndHex0,psupHex0,aQuadHex0,
      aHex0.data(), aHex0.shape()[0],
      aXYZ0.data(), aXYZ0.shape()[0]);
  py::array_t<double> npXYZ1({(int)aXYZ1.size()/3,3}, aXYZ1.data());
  py::array_t<unsigned int> npHex1({(int)aHex1.size()/8,8}, aHex1.data());
  return std::make_tuple(npXYZ1,npHex1);
}

// -------------------------------------------------------------

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

std::tuple<
    py::array_t<unsigned int>,
    py::array_t<unsigned int>
    >
PyJArray_MeshPsup(
    const py::array_t<unsigned int>& elm,
    int npoint)
{
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_PSuP_MeshElem(psup_ind, psup,
                             elm.data(), elm.shape()[0], elm.shape()[1],
                             npoint);
  py::array_t<unsigned int> np_psup_ind((pybind11::size_t)psup_ind.size(), psup_ind.data());
  py::array_t<unsigned int> np_psup((pybind11::size_t)psup.size(), psup.data());
  return std::make_tuple(np_psup_ind, np_psup);
}

void
PyJArray_Sort(
    const py::array_t<unsigned int>& psup_ind,
    py::array_t<unsigned int>& psup)
{
  assert( dfm2::CheckJArray(psup_ind,psup) );
  dfm2::JArray_Sort(psup_ind.data(), psup_ind.shape(0)-1, psup.mutable_data() );
}

std::tuple<
    py::array_t<unsigned int>,
    py::array_t<unsigned int>
    >
PyJArray_Extend(
    const py::array_t<unsigned int>& psup_ind,
    const py::array_t<unsigned int>& psup)
{
  assert( dfm2::CheckJArray(psup_ind,psup) );
  std::vector<unsigned int> psup_ind1, psup1;
  dfm2::JArray_Extend(
      psup_ind1, psup1,
      psup_ind.data(), psup_ind.shape()[0], psup.data() );
  py::array_t<unsigned int> np_psup_ind1((pybind11::size_t)psup_ind1.size(), psup_ind1.data());
  py::array_t<unsigned int> np_psup1((pybind11::size_t)psup1.size(), psup1.data());
  return std::make_tuple(np_psup_ind1, np_psup1);
}

std::tuple<
    py::array_t<unsigned int>,
    py::array_t<unsigned int>
    >
PyJArray_AddDiagonal(
    const py::array_t<unsigned int>& psup_ind0,
    const py::array_t<unsigned int>& psup0)
{
  std::vector<unsigned int> psup_ind, psup;
  dfm2::JArray_AddDiagonal(psup_ind,psup,
                           psup_ind0.data(),psup_ind0.shape()[0], psup0.data(),psup0.shape()[0]);
  py::array_t<unsigned int> np_psup_ind((pybind11::size_t)psup_ind.size(), psup_ind.data());
  py::array_t<unsigned int> np_psup((pybind11::size_t)psup.size(), psup.data());
  return std::make_tuple(np_psup_ind, np_psup);
}

py::array_t<unsigned int>
PyElemQuad_DihedralTri(
    const py::array_t<unsigned int>& aTri,
    int np)
{
  assert( dfm2::CheckNumpyArray2D(aTri, -1, 3) );
  const unsigned int nTri = aTri.shape()[0];
  std::vector<unsigned int> aQuad;
  dfm2::ElemQuad_DihedralTri(aQuad, aTri.data(), nTri, np);
  py::array_t<unsigned int> npQuad({(int)aQuad.size()/4,4}, aQuad.data());
  return npQuad;
}

std::tuple<double,double>
PyQuality_MeshTri2D(
    const py::array_t<double>& np_xy,
    const py::array_t<unsigned int>& np_tri)
{
  assert( dfm2::CheckNumpyArray2D(np_xy, -1, 2) );
  assert( dfm2::CheckNumpyArray2D(np_tri, -1, 3) );
  double max_aspect;
  double min_area;
  dfm2::Quality_MeshTri2D(max_aspect, min_area,
      np_xy.data(),
      np_tri.data(), np_tri.shape()[0]);
   return std::make_tuple(max_aspect, min_area);
}

py::array_t<double>
PyMassPoint_MeshTri(
    py::array_t<double>& npXY,
    py::array_t<unsigned int>& npTri)
{
  assert( dfm2::CheckNumpyArray2D(npTri, -1, 3) );
  assert( npXY.ndim() == 2 );
  const unsigned int np = npXY.shape()[0];
  py::array_t<double> npMass(np);
  if( npXY.shape()[1] == 2 ) {
    dfm2::MassPoint_Tri2D(npMass.mutable_data(), 1.0,
                          npXY.data(), npXY.shape()[0],
                          npTri.data(), npTri.shape()[0]);
  }
  else if( npXY.shape()[1] == 3 ) {
    dfm2::MassPoint_Tri3D(npMass.mutable_data(), 1.0,
                          npXY.data(), npXY.shape()[0],
                          npTri.data(), npTri.shape()[0]);
  }
  return npMass;
}

void
PyNormalVtx_Mesh(
    py::array_t<double>& npNrm,
    const py::array_t<double>& pos,
    const py::array_t<unsigned int>& elm,
    const dfm2::MESHELEM_TYPE type)
{
  assert( dfm2::CheckNumpyArray2D(pos, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(elm, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npNrm, -1, 3) );
  assert( npNrm.shape()[0] == pos.shape()[0] );
  if( type == dfm2::MESHELEM_TRI ){
    delfem2::Normal_MeshTri3D(npNrm.mutable_data(),
                              pos.data(),pos.shape()[0],
                              elm.data(),elm.shape()[0]);
  }
}

py::array_t<unsigned int>
PyEdge_Mesh(
    const py::array_t<double>& pos,
    const py::array_t<unsigned int>& elm,
    const dfm2::MESHELEM_TYPE type)
{
  assert( dfm2::CheckNumpyArray2D(pos, -1, pos.shape()[1]) );
  assert( dfm2::CheckNumpyArray2D(elm, -1, nNodeElem(type)) );
  std::vector<unsigned int> elsup_ind, elsup;
  dfm2::JArray_ElSuP_MeshElem(elsup_ind, elsup,
                              elm.data(), elm.shape()[0], elm.shape()[1], pos.shape()[0]);
  std::vector<unsigned int> edge_ind, edge;
  JArrayEdge_MeshElem(edge_ind, edge,
                      elm.data(), type, elsup_ind, elsup, false);
  std::vector<unsigned int> aLine;
  dfm2::MeshLine_JArrayEdge(aLine,
                            edge_ind, edge);
  py::array_t<unsigned int> npLine({(int)aLine.size()/2,2}, aLine.data());
  return npLine;
}

// --------------------------------------

void init_mshtopoio(py::module &m){
  py::enum_<dfm2::MESHELEM_TYPE>(m, "MESH_ELEM_TYPE")
  .value("TRI",     dfm2::MESHELEM_TYPE::MESHELEM_TRI)
  .value("QUAD",    dfm2::MESHELEM_TYPE::MESHELEM_QUAD)
  .value("TET",     dfm2::MESHELEM_TYPE::MESHELEM_TET)
  .value("PYRAMID", dfm2::MESHELEM_TYPE::MESHELEM_PYRAMID)
  .value("WEDGE",   dfm2::MESHELEM_TYPE::MESHELEM_WEDGE)
  .value("HEX",     dfm2::MESHELEM_TYPE::MESHELEM_HEX)
  .value("LINE",    dfm2::MESHELEM_TYPE::MESHELEM_LINE)
  .export_values();
  
  py::class_<CMeshMultiElem>(m,"CppMeshMultiElem")
  .def(py::init<>())
  .def("minmax_xyz",  &CMeshMultiElem::AABB3_MinMax)
  //
  .def("read_obj",    &CMeshMultiElem::ReadObj)
  .def("scaleXYZ",    &CMeshMultiElem::ScaleXYZ)
  .def("translateXYZ",&CMeshMultiElem::TranslateXYZ);

  m.def("num_node_elem", &dfm2::nNodeElem);

  m.def("cppNormalVtx_Mesh",      &PyNormalVtx_Mesh);
  m.def("cppEdge_Mesh",           &PyEdge_Mesh);
  
  // io read
  m.def("meshtri3d_read_ply",     &PyMeshTri3D_ReadPly,     py::return_value_policy::move);
  m.def("meshtri3d_read_obj",     &PyMeshTri3D_ReadObj,     py::return_value_policy::move);
  m.def("meshtri3d_read_nastran", &PyMeshTri3D_ReadNastran, py::return_value_policy::move);

  // io write
  m.def("meshtri3d_write_obj",    &PyMeshTri3D_WriteObj);
  
 // subdiv
  m.def("meshquad3d_subdiv",      &PyMeshQuad3D_Subviv,     py::return_value_policy::move);
  m.def("meshhex3d_subdiv",       &PyMeshHex3D_Subviv,      py::return_value_policy::move);
  
  m.def("setTopology_ExtrudeTri2Tet", &PySetTopology_ExtrudeTri2Tet);
  
  // jarray
  m.def("cppJArray_MeshPsup",    &PyJArray_MeshPsup,    py::return_value_policy::move);
  m.def("cppJarray_AddDiagonal", &PyJArray_AddDiagonal, py::return_value_policy::move);
  m.def("cppJArray_Sort",        &PyJArray_Sort,        py::return_value_policy::move);
  m.def("cppJArray_Extend",      &PyJArray_Extend,      py::return_value_policy::move);

  m.def("elemQuad_dihedralTri",&PyElemQuad_DihedralTri);
  m.def("quality_meshTri2D",   &PyQuality_MeshTri2D);
  m.def("cppMassPoint_MeshTri",&PyMassPoint_MeshTri);
}
