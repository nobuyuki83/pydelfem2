/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "delfem2/tinygltf/io_gltf.h"
#include "delfem2/rig_geo3.h"
#include "../py_funcs.h"

#include "tinygltf/tiny_gltf.h"
#include "stb_image.h" // stb is already compiled in io_gltf.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <deque>


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

class CSkeleton{
public:
  void SetTranslation(int ib, const std::vector<double>& aT){
    assert(aT.size()==3);
    aRigBone[ib].SetTranslation(aT[0], aT[1], aT[2]);
    UpdateBoneRotTrans(aRigBone);
  }
  void SetRotationBryant(int ib, const std::vector<double>& aRB){
    assert(aRB.size()==3);
    double dq[4]; dfm2::Quat_Bryant(dq,aRB[0],aRB[1],aRB[2]);
    double q0[4]; dfm2::Copy_Quat(q0, aRigBone[ib].quatRelativeRot);
    dfm2::QuatQuat(aRigBone[ib].quatRelativeRot,q0,dq);
//    aRigBone[ib].SetRotationBryant(aRB[0], aRB[1], aRB[2]);
    UpdateBoneRotTrans(aRigBone);
  }
public:
  std::vector<dfm2::CRigBone> aRigBone;
};

CSkeleton
PySkeleton_Gltf(
    const dfm2::CGLTF& gltf,
    int iskin)
{
  CSkeleton bones;
  gltf.GetBone(bones.aRigBone,
               iskin);
  return bones;
}

CSkeleton
PyGetSkeleton_BoneTreePos(
  const py::array_t<unsigned int>& npIdBoneParent,
  const py::array_t<double>& npJ)
{
  assert( npIdBoneParent.ndim() == 1 && npIdBoneParent.strides(0) == sizeof(unsigned int) );
  assert( dfm2::CheckNumpyArray2D(npJ, -1, 3) );
  const unsigned int nb = npIdBoneParent.shape()[0];
  assert( npJ.shape()[0] == nb );
  CSkeleton skeleton;
  dfm2::InitBones_JointPosition(skeleton.aRigBone,
      nb, npIdBoneParent.data(), npJ.data());
  return skeleton;
}

void PyUpdateRigSkin(
    py::array_t<double>& npXYZ,
    const py::array_t<double>& npXYZ0,
    const CSkeleton& skeleton,
    const py::array_t<double>& npRigWeight,
    const py::array_t<unsigned int>& npRigJoint)
{
  assert( dfm2::CheckNumpyArray2D(npXYZ, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npXYZ0, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npRigWeight, -1, 4) );
  assert( dfm2::CheckNumpyArray2D(npRigJoint, -1, 4) );
  assert( npXYZ.shape()[0] == npXYZ0.shape()[0] );
  assert( npXYZ.shape()[0] == npRigWeight.shape()[0] );
  assert( npXYZ.shape()[0] == npRigJoint.shape()[0] );
  dfm2::Skinning_LBS_LocalWeight(npXYZ.mutable_data(),
                                 npXYZ0.data(), npXYZ0.shape()[0],
                                 skeleton.aRigBone,
                                 npRigWeight.data(),
                                 npRigJoint.data());
}

std::tuple< py::array_t<double>,py::array_t<unsigned int> >
PySparsifyMatrixRow(
    const py::array_t<double>& A)
{
  std::vector<double> aSparseA;
  std::vector<unsigned int> aSparseI;
  const int nrow = A.shape()[0];
  dfm2::SparsifyMatrixRow(aSparseA,aSparseI,A.data(),nrow,A.shape()[1],1.0e-4);
  const int nelm_par_row = aSparseA.size()/nrow;
  py::array_t<unsigned int> npSparseI({nrow,nelm_par_row}, aSparseI.data());
  py::array_t<double> npSparseA({nrow,nelm_par_row}, aSparseA.data());
  return std::make_tuple(npSparseA,npSparseI);
}


// Rigging related ends here
// -----------------------------------------

void init_rigging(py::module &m){

  py::class_<dfm2::CGLTF>(m,"CppGLTF")
      .def(py::init<>())
      .def("read", &dfm2::CGLTF::Read)
      .def("print", &dfm2::CGLTF::Print);

  // I defined the array as a class because pybind11 cannot use call by reference for list.
  py::class_<CSkeleton>(m,"CppSkeleton")
      .def("set_translation", &CSkeleton::SetTranslation)
      .def("set_rotation_bryant", &CSkeleton::SetRotationBryant)
      .def(py::init<>());

  m.def("cppGetMeshInfoGltf", &PyGLTF_GetMeshInfo);
  m.def("cppGetSkeleton_Gltf", &PySkeleton_Gltf);
  m.def("cppUpdateRigSkin", &PyUpdateRigSkin);
  m.def("cppUpdateBoneRotTransl", &dfm2::UpdateBoneRotTrans);
  m.def("cppGetSkeleton_BoneTreePos",&PyGetSkeleton_BoneTreePos);
  m.def("cppSparsifyMatrixRow",&PySparsifyMatrixRow);


}
