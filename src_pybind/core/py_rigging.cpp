/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <vector>
#include <deque>

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/rig_geo3.h"

#include "../py_funcs.h"
#include "tinygltf/tiny_gltf.h"

#include "stb_image.h" // stb is already compiled in io_gltf.cpp

namespace py = pybind11;
namespace dfm2 = delfem2;

// -----------------------------------------------
// Rigging related from here

std::tuple<py::array_t<double>, py::array_t<unsigned int>, py::array_t<double>, py::array_t<unsigned int>>
PyGLTF_GetMeshInfo(
    const dfm2::CGLTF& gltf,
    int imesh,
    int iprimitive)
{
  std::vector<double> aXYZ0;
  std::vector<unsigned int> aTri;
  std::vector<double> aRigWeight;
  std::vector<unsigned int> aRigJoint;
  gltf.GetMeshInfo(aXYZ0,aTri,aRigWeight,aRigJoint,
                   imesh, iprimitive);
  const unsigned int np = aXYZ0.size()/3;
  assert( aRigWeight.size() == np*4 );
  assert( aRigJoint.size() == np*4 );
  py::array_t<double> npXYZ0({(int)np,3}, aXYZ0.data());
  py::array_t<unsigned int> npTri({(int)aTri.size()/3,3}, aTri.data());
  py::array_t<double> npRW({(int)np,4}, aRigWeight.data());
  py::array_t<unsigned int> npRJ({(int)np,4}, aRigJoint.data());
  return std::make_tuple(npXYZ0,npTri,npRW,npRJ);
}

class CBoneArray{
public:
  void SetTranslation(int ib, const std::vector<double>& aT){
    assert(aT.size()==3);
    aRigBone[ib].SetTranslation(aT[0], aT[1], aT[2]);
    UpdateBoneRotTrans(aRigBone);
  }
  void SetRotationBryant(int ib, const std::vector<double>& aRB){
    assert(aRB.size()==3);
    aRigBone[ib].SetRotationBryant(aRB[0], aRB[1], aRB[2]);
    UpdateBoneRotTrans(aRigBone);
  }
public:
  std::vector<dfm2::CRigBone> aRigBone;
};

CBoneArray
PyGLTF_GetBones
(const dfm2::CGLTF& gltf,
 int iskin)
{
  CBoneArray BA;
  gltf.GetBone(BA.aRigBone,
               iskin);
  return BA;
}

void PyUpdateRigSkin(
    py::array_t<double>& npXYZ,
    const py::array_t<double>& npXYZ0,
    const py::array_t<unsigned int>& npTri,
    const CBoneArray& BA,
    const py::array_t<double>& npRigWeight,
    const py::array_t<unsigned int>& npRigJoint)
{
  assert( dfm2::CheckNumpyArray2D(npXYZ, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npXYZ0, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npTri, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npRigWeight, -1, 4) );
  assert( dfm2::CheckNumpyArray2D(npRigJoint, -1, 4) );
  assert( npXYZ.shape()[0] == npXYZ0.shape()[0] );
  assert( npXYZ.shape()[0] == npRigWeight.shape()[0] );
  assert( npXYZ.shape()[0] == npRigJoint.shape()[0] );
  dfm2::Skinning_LBS_LocalWeight(npXYZ.mutable_data(),
                                 npXYZ0.data(), npXYZ0.shape()[0],
                                 npTri.data(), npTri.shape()[0],
                                 BA.aRigBone,
                                 npRigWeight.data(),
                                 npRigJoint.data());
}

// Rigging related ends here
// -----------------------------------------

void init_rigging(py::module &m){

  py::class_<dfm2::CGLTF>(m,"CppGLTF")
      .def(py::init<>())
      .def("read", &dfm2::CGLTF::Read)
      .def("print", &dfm2::CGLTF::Print);

  py::class_<CBoneArray>(m,"CppBoneArray")
      .def("set_translation", &CBoneArray::SetTranslation)
      .def("set_rotation_bryant", &CBoneArray::SetRotationBryant)
      .def(py::init<>());

  m.def("cppGetMeshInfoGltf", &PyGLTF_GetMeshInfo);
  m.def("cppGetBonesGltf", &PyGLTF_GetBones);
  m.def("cppUpdateRigSkin", &PyUpdateRigSkin);
  m.def("cppUpdateBoneRotTransl", &dfm2::UpdateBoneRotTrans);

}
