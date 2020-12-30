/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "delfem2/mshtopoio.h" // for GeodesicPolygon
#include "delfem2/mshprimitive.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// ----------------------------------------------------------

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshQuad2D_Grid(
    int mx,
    int my)
{
  std::vector<double> aXY;
  std::vector<unsigned int> aQuad;
  dfm2::MeshQuad2D_Grid(aXY, aQuad,
                        mx-1, my-1);
  py::array_t<double> npXY({(int)aXY.size()/2,2}, aXY.data());
  py::array_t<unsigned int> npQuad({(int)aQuad.size()/4,4}, aQuad.data());
  return std::make_tuple(npXY,npQuad);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshTri3D_Cylinder(
    double r,
    double l,
    int nr,
    int nl)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_CylinderClosed(aXYZ, aTri, r, l, nr, nl);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}

std::tuple<
    py::array_t<double>,
    py::array_t<unsigned int>
    >
PyMeshTri3D_Cube(int n)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Cube(aXYZ, aTri, n);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}

std::tuple< py::array_t<double>, py::array_t<unsigned int> >
PyMeshTri3D_Sphere(double r, int nla, int nlo)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Sphere(aXYZ,aTri, r,nla,nlo);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}

std::tuple< py::array_t<double>, py::array_t<unsigned int> >
PyMeshTri3D_GeoPoly()
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  MeshTri3D_GeodesicPolyhedron(aXYZ,aTri);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}


std::tuple< py::array_t<double>, py::array_t<unsigned int> >
PyMeshTri3D_Icosahedron()
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3D_Icosahedron(aXYZ,aTri);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}


std::tuple< py::array_t<double>, py::array_t<unsigned int> >
PyMeshTri3D_Torus(double r0, double r1)
{
  std::vector<double> aXYZ;
  std::vector<unsigned int> aTri;
  dfm2::MeshTri3_Torus(aXYZ,aTri,
                       r0, r1, 18, 32);
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npXYZ({(int)aXYZ.size()/3,3}, aXYZ.data());
  return std::make_tuple(npXYZ,npTri);
}

// above: primitive
// --------------------------------------

void init_primitive(py::module &m){
  // primitive
  m.def("cppMeshQuad2D_Grid",       &PyMeshQuad2D_Grid,       py::return_value_policy::move);
  m.def("cppMeshTri3D_Torus",       &PyMeshTri3D_Torus,       py::return_value_policy::move);
  m.def("cppMeshTri3D_Cylinder",    &PyMeshTri3D_Cylinder,    py::return_value_policy::move);
  m.def("cppMeshTri3D_Sphere",      &PyMeshTri3D_Sphere,      py::return_value_policy::move);
  m.def("cppMeshTri3D_Cube",        &PyMeshTri3D_Cube,        py::return_value_policy::move);
  m.def("cppMeshTri3D_GeoPoly",     &PyMeshTri3D_GeoPoly,     py::return_value_policy::move);
  m.def("cppMeshTri3D_Icosahedron", &PyMeshTri3D_Icosahedron, py::return_value_policy::move);
}
