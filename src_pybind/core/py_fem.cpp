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
#include "delfem2/mshuni.h"
#include "delfem2/mshmisc.h"
#include "delfem2/sdf.h"

#include "delfem2/dtri2_v2dtri.h"

#include "delfem2/fempoisson.h"
#include "delfem2/femsolidlinear.h"
#include "delfem2/femstokes.h"
#include "delfem2/femnavierstokes.h"
#include "delfem2/femmips_geo3.h"
#include "delfem2/femmitc3.h"
#include "delfem2/femcloth.h"
#include "delfem2/pbd_geo3dtri23.h"

#include "delfem2/jagarray.h"

namespace py = pybind11;
namespace dfm2 = delfem2;

// ------------------------------------------------------------

void PyMergeLinSys_Poission
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double alpha, double source,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aElm,
     dfm2::MESHELEM_TYPE elem_type,
     const py::array_t<double>& aVal)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == dfm2::MESHELEM_TRI ){
      dfm2::MergeLinSys_Poission_MeshTri2D(mss, (double*)buff_vecb.ptr,
                                           alpha, source,
                                           aXY.data(), aXY.shape()[0],
                                           aElm.data(), aElm.shape()[0],
                                           aVal.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == dfm2::MESHELEM_TET ){
      dfm2::MergeLinSys_Poission_MeshTet3D(mss, (double*)buff_vecb.ptr,
                                           alpha, source,
                                           aXY.data(), aXY.shape()[0],
                                           aElm.data(), aElm.shape()[0],
                                           aVal.data());
    }
  }
}


void PyMergeLinSys_Diffuse
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double alpha, double rho, double source,
     double dt_timestep, double gamma_newmark,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aElm,
     dfm2::MESHELEM_TYPE elem_type,
     const py::array_t<double>& aVal,
     const py::array_t<double>& aVelo)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == dfm2::MESHELEM_TRI ){
      dfm2::MergeLinSys_Diffusion_MeshTri2D(mss, (double*)buff_vecb.ptr,
                                            alpha, rho, source,
                                            dt_timestep, gamma_newmark,
                                            aXY.data(), aXY.shape()[0],
                                            aElm.data(), aElm.shape()[0],
                                            aVal.data(), aVelo.data());
    }
  }
  else if( aXY.shape()[1] == 3 ){
    if( elem_type == dfm2::MESHELEM_TET ){
      dfm2::MergeLinSys_Diffusion_MeshTet3D(mss, (double*)buff_vecb.ptr,
                                            alpha, rho, source,
                                            dt_timestep, gamma_newmark,
                                            aXY.data(), aXY.shape()[0],
                                            aElm.data(), aElm.shape()[0],
                                            aVal.data(), aVelo.data());
    }
  }
}




void PyMergeLinSys_LinearSolidStatic
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double myu, double lambda, double rho,
     std::vector<double>& gravity,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aElm,
     dfm2::MESHELEM_TYPE elem_type,
     const py::array_t<double>& aVal)
{
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  auto buff_vecb = vec_b.request();
  if( aXY.shape()[1] == 2 ){
    if( elem_type == dfm2::MESHELEM_TRI ){
      dfm2::MergeLinSys_SolidLinear_Static_MeshTri2D(
          mss, (double*)buff_vecb.ptr,
          myu,lambda,rho,
          gravity[0], gravity[1],
          aXY.data(), aXY.shape()[0],
          aElm.data(), aElm.shape()[0],
          aVal.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == dfm2::MESHELEM_TET ){
      dfm2::MergeLinSys_SolidLinear_Static_MeshTet3D(
          mss, (double*)buff_vecb.ptr,
          myu,lambda,rho,gravity.data(),
          aXY.data(), aXY.shape()[0],
          aElm.data(), aElm.shape()[0],
          aVal.data());
    }
  }
}


void PyMergeLinSys_LinearSolidDynamic
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double myu, double lambda, double rho,
     std::vector<double>& gravity,
     double dt_timestep, double gamma_newmark, double beta_newmark,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aElm,
     dfm2::MESHELEM_TYPE elem_type,
     const py::array_t<double>& aVal,
     const py::array_t<double>& aVelo,
     const py::array_t<double>& aAcc)
{
  auto buff_vecb = vec_b.request();
  assert( aXY.shape()[1] == 2 || aXY.shape()[1] == 3 );
  assert( nNodeElem(elem_type) == aElm.shape()[1] );
  if( aXY.shape()[1] == 2 ){
    if( elem_type == dfm2::MESHELEM_TRI ){
      dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTri2D(
          mss,(double*)buff_vecb.ptr,
          myu,lambda,rho,gravity[0],gravity[1],
          dt_timestep,gamma_newmark,beta_newmark,
          aXY.data(), aXY.shape()[0],
          aElm.data(), aElm.shape()[0],
          aVal.data(),aVelo.data(),aAcc.data());
    }
  }
  if( aXY.shape()[1] == 3 ){
    if( elem_type == dfm2::MESHELEM_TET ){
      dfm2::MergeLinSys_SolidLinear_NewmarkBeta_MeshTet3D(
          mss,(double*)buff_vecb.ptr,
          myu,lambda,rho,gravity.data(),
          dt_timestep,gamma_newmark,beta_newmark,
          aXY.data(), aXY.shape()[0],
          aElm.data(), aElm.shape()[0],
          aVal.data(),aVelo.data(),aAcc.data());
    }
  }
}


void PyMergeLinSys_StorksStatic2D
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double myu, double g_x, double g_y,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aTri,
     const py::array_t<double>& aVal)
{
  auto buff_vecb = vec_b.request();
  dfm2::MergeLinSys_StokesStatic2D(
      mss,(double*)buff_vecb.ptr,
      myu,g_x,g_y,
      aXY.data(), aXY.shape()[0],
      aTri.data(), aTri.shape()[0],
      aVal.data());
}

void PyMergeLinSys_StorksDynamic2D
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double myu, double rho, double g_x, double g_y,
     double dt_timestep, double gamma_newmark,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aTri,
     const py::array_t<double>& aVal,
     const py::array_t<double>& aVelo)
{
  auto buff_vecb = vec_b.request();
  dfm2::MergeLinSys_StokesDynamic2D(
      mss,(double*)buff_vecb.ptr,
      myu,rho,g_x,g_y,
      dt_timestep, gamma_newmark,
      aXY.data(), aXY.shape()[0],
      aTri.data(), aTri.shape()[0],
      aVal.data(),aVelo.data());
}

void PyMergeLinSys_NavierStorks2D
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double myu, double rho, double g_x, double g_y,
     double dt_timestep, double gamma_newmark,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aTri,
     const py::array_t<double>& aVal,
     const py::array_t<double>& aVelo)
{
  auto buff_vecb = vec_b.request();
  dfm2::MergeLinSys_NavierStokes2D(
      mss,(double*)buff_vecb.ptr,
      myu,rho,g_x,g_y,
      dt_timestep, gamma_newmark,
      aXY.data(), aXY.shape()[0],
      aTri.data(), aTri.shape()[0],
      aVal.data(), aVelo.data());
}

void PyMergeLinSys_ShellMitc3Static
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double thickness, double lambda, double myu, double rho, double g_z,
     const py::array_t<double>& aXY,
     const py::array_t<unsigned int>& aTri,
     const py::array_t<double>& aDisp)
{
  auto buff_vecb = vec_b.request();
  dfm2::MergeLinSys_ShellStaticPlateBendingMITC3_MeshTri2D(
      mss,(double*)buff_vecb.ptr,
      thickness, lambda, myu, rho, g_z,
      aXY.data(), aXY.shape()[0],
      aTri.data(), aTri.shape()[0],
      aDisp.data());
}

double PyMergeLinSys_Cloth
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double lambda, double myu, double stiff_bend,
     const py::array_t<double>& aPosIni,
     const py::array_t<unsigned int>& aTri,
     const py::array_t<unsigned int>& aQuad,
     const py::array_t<double>& aXYZ)
{
  auto buff_vecb = vec_b.request();
  double W = dfm2::MergeLinSys_Cloth(
      mss,(double*)buff_vecb.ptr,
      lambda, myu, stiff_bend,
      aPosIni.data(), aPosIni.shape()[0], aPosIni.shape()[1],
      aTri.data(), aTri.shape()[0],
      aQuad.data(), aQuad.shape()[0],
      aXYZ.data());
  return W;
}


double PyMergeLinSys_Contact
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
        ////
     double stiff_contact,
     double contact_clearance,
     const std::vector<const dfm2::CSDF3*>& apSDF,
     const py::array_t<double>& aXYZ)
{
  if( apSDF.empty() ) return 0;
  class CMyInput : public dfm2::CInput_Contact
  {
  public:
    explicit CMyInput(const std::vector<const dfm2::CSDF3*>& apSDF){ this->apSDF = apSDF; }
    double penetrationNormal(double& nx, double& ny, double& nz,
                             double px, double py, double pz) const override
    {
      double n[3];
      double max_pd = apSDF[0]->Projection(n,
                                           px,py,pz);
      /*
      for(unsigned int ipct=1;ipct<apSDF.size();ipct++){
        double dist0,n0[3];
        dist0 = apSDF[ipct]->Projection(px,py,pz, n0);
        if( dist0 < max_pd ) continue;
        max_pd = dist0;
        n[0] = n0[0];
        n[1] = n0[1];
        n[2] = n0[2];
      }
       */
      nx = -n[0];
      ny = -n[1];
      nz = -n[2];
      return max_pd;
    }
  public:
    std::vector<const dfm2::CSDF3*> apSDF;
  } input(apSDF);
  auto buff_vecb = vec_b.request();
  double W = dfm2::MergeLinSys_Contact(
      mss, (double*)buff_vecb.ptr,
      stiff_contact,contact_clearance,
      input,
      aXYZ.data(), aXYZ.shape()[0]);
  return W;
}

double PyMergeLinSys_MassPoint
    (dfm2::CMatrixSparse<double>& mss,
     py::array_t<double>& vec_b,
     double mass_point,
     double dt,
     const std::vector<double>& gravity,
     const py::array_t<double>& aXYZ,
     const py::array_t<double>& aUVW)
{
  double* pB = (double*)(vec_b.request().ptr);
  double W = 0.0;
  const int np = aXYZ.shape()[0];
  assert(aUVW.shape()[0] == np);
  for(int ip=0;ip<np;ip++){
    const double c[3] = {aXYZ.at(ip,0),aXYZ.at(ip,1),aXYZ.at(ip,2)};
    W -= mass_point*( c[0]*gravity[0] + c[1]*gravity[1] + c[2]*gravity[2] );
    pB[ip*3+0] -= mass_point*gravity[0];
    pB[ip*3+1] -= mass_point*gravity[1];
    pB[ip*3+2] -= mass_point*gravity[2];
  }
  const int ndof = aXYZ.size();
  const double* pUVW = aUVW.data();
  for(int i=0;i<ndof;i++){
    pB[i] = -pB[i] + mass_point*pUVW[i]/dt;
  }
  for(int ip=0;ip<np;ip++){
    mss.valDia[ip*9+0*3+0] += mass_point / (dt*dt);
    mss.valDia[ip*9+1*3+1] += mass_point / (dt*dt);
    mss.valDia[ip*9+2*3+2] += mass_point / (dt*dt);
  }
  return W;
}

void PyPBD_ConstProj_Rigid2D
    (py::array_t<double>& npXYt,
     double stiffness,
     const py::array_t<unsigned int>& npClstrInd,
     const py::array_t<unsigned int>& npClstr,
     const py::array_t<double>& npXY)
{
  dfm2::PBD_ConstProj_Rigid2D((double*)(npXYt.request().ptr),
                              stiffness,
                              npClstrInd.data(), npClstrInd.size(),
                              npClstr.data(),    npClstr.size(),
                              npXY.data(),       npXY.shape()[0]);
}

void PyConstProj_Rigid3D(
    py::array_t<double>& npXYZt,
    double stiffness,
    const py::array_t<int>& npClstrInd,
    const py::array_t<int>& npClstr,
    const py::array_t<double>& npXYZ)
{
  assert( dfm2::CheckNumpyArray2D(npXYZt,-1,3) );
  dfm2::PBD_ConstProj_Rigid3D(npXYZt.mutable_data(),
      stiffness,
      npClstrInd.data(), npClstrInd.size(),
      npClstr.data(),    npClstr.size(),
      npXYZ.data(),      npXYZ.shape()[0]);
}

void PyPBD_ConstProj_ClothStretch
    (py::array_t<double>& npXYZt,
     const dfm2::CMeshDynTri2D& mesh)
{
  assert( dfm2::CheckNumpyArray2D(npXYZt,-1,3) );
  const std::vector<dfm2::CDynTri>& aETri = mesh.aETri;
  const std::vector<dfm2::CVec2d>& aVec2 = mesh.aVec2;
  PBD_TriStrain(npXYZt.mutable_data(),
      npXYZt.shape()[0], aETri, aVec2);
}

void PyPBD_ConstProj_ClothBend
    (py::array_t<double>& npXYZt,
     const dfm2::CMeshDynTri2D& mesh)
{
  assert( dfm2::CheckNumpyArray2D(npXYZt,-1,3) );
  const std::vector<dfm2::CDynTri>& aETri = mesh.aETri;
  const std::vector<dfm2::CVec2d>& aVec2 = mesh.aVec2;
  PBD_Bend(npXYZt.mutable_data(),
      npXYZt.shape()[0],
      aETri, aVec2,
      1.0);
}


void PyPBD_ConstProj_Seam
    (py::array_t<double>& npXYZt,
     const py::array_t<unsigned int>& npLine)
{
  assert( dfm2::CheckNumpyArray2D(npXYZt, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(npLine, -1, 2) );
  const unsigned int nline = npLine.shape()[0];
  dfm2::PBD_Seam(npXYZt.mutable_data(),
      npXYZt.shape()[0],
      npLine.data(), nline);
}

void PyPBD_ConstProj_Contact
    (py::array_t<double>& npXYZt,
     const dfm2::CSDF3& sdf)
{
  assert( dfm2::CheckNumpyArray2D(npXYZt, -1, 3) );
  double* aXYZt = npXYZt.mutable_data();
  const unsigned int np = npXYZt.shape()[0];
  for(unsigned int ip=0;ip<np;++ip){
    double n[3];
    double dist = sdf.Projection(n,
                                 aXYZt[ip*3+0], aXYZt[ip*3+1], aXYZt[ip*3+2]);
    if( dist > 0 ){
      aXYZt[ip*3+0] += dist*n[0];
      aXYZt[ip*3+1] += dist*n[1];
      aXYZt[ip*3+2] += dist*n[2];
    }
  }
//  assert( npLine.ndim() == 2 );
//  assert( npLine.shape()[1] == 2 );
//  const unsigned int nline = npLine.shape()[0];
}

void PyPointFixBC
    (py::array_t<double>& aTmp,
     const py::array_t<int>& aBC,
     const py::array_t<double>& npXY1)
{
  assert( aTmp.ndim() == 2 );
  assert( npXY1.ndim() == 2 );
  assert( aTmp.shape()[1] == npXY1.shape()[1] );
  const int np = aTmp.shape()[0];
  double* ptr = aTmp.mutable_data();
  if( npXY1.shape()[1] == 2 ){
    for(int ip=0;ip<np;++ip){
      if( aBC.at(ip) == 0 ){ continue; }
      ptr[ip*2+0] = npXY1.at(ip,0);
      ptr[ip*2+1] = npXY1.at(ip,1);
    }
  }
  if( npXY1.shape()[1] == 3 ){
    for(int ip=0;ip<np;++ip){
      if( aBC.at(ip) == 0 ){ continue; }
      ptr[ip*3+0] = npXY1.at(ip,0);
      ptr[ip*3+1] = npXY1.at(ip,1);
      ptr[ip*3+2] = npXY1.at(ip,2);
    }
  }
}


void PyMassPointMesh
    (py::array_t<double>& mass_point,
     double rho,
     const py::array_t<double>& np_pos,
     const py::array_t<unsigned int>& np_elm,
     dfm2::MESHELEM_TYPE elem_type)
{
  assert( mass_point.ndim() == 1 );
  assert( np_pos.ndim() == 2 );
  assert( np_elm.ndim() == 2 );
  assert( mass_point.shape()[0] == np_pos.shape()[0] );
  assert( dfm2::CheckNumpyArray2D(np_elm, -1, nNodeElem(elem_type)) );
  if( elem_type ==  dfm2::MESHELEM_TET ){
    assert( dfm2::CheckNumpyArray2D(np_pos, -1, 3) );
    dfm2::MassPoint_Tet3D(mass_point.mutable_data(),
                          rho,
                          np_pos.data(), np_pos.shape()[0],
                          np_elm.data(), np_elm.shape()[0]);
  }
  else if( elem_type ==  dfm2::MESHELEM_TRI ){
    if( np_pos.shape()[1] == 2 ){ // two dimensional
      assert( dfm2::CheckNumpyArray2D(np_pos, -1, 2) );
      dfm2::MassPoint_Tri2D(mass_point.mutable_data(),
                            rho,
                            np_pos.data(), np_pos.shape()[0],
                            np_elm.data(), np_elm.shape()[0]);
    }
    else{
      assert(0);
    }
  }
  else{
    // TODO: implemnet mass lumped for other types of meshes
    assert(0);
  }
}

void PyMassLumped_ShellPlateBendingMitc3
    (py::array_t<double>& mass_lumped,
     double rho, double thick,
     const py::array_t<double>& np_pos,
     const py::array_t<unsigned int>& np_elm)
{
  assert( dfm2::CheckNumpyArray2D(mass_lumped, -1, 3) );
  assert( dfm2::CheckNumpyArray2D(np_pos, -1, 2) );
  assert( dfm2::CheckNumpyArray2D(np_elm, -1, 3) );
  assert( mass_lumped.shape()[0] == np_pos.shape()[0] );
  dfm2::MassLumped_ShellPlateBendingMITC3(mass_lumped.mutable_data(),
      rho, thick,
      np_pos.data(), np_pos.shape()[0],
      np_elm.data(), np_elm.shape()[0]);
}

void init_fem(py::module &m)
{
  m.def("cppMassPoint_Mesh",                    &PyMassPointMesh);
  m.def("cppMassLumped_ShellPlateBendingMitc3", &PyMassLumped_ShellPlateBendingMitc3);

  m.def("cppFEM_Merge_PointMass",         &PyMergeLinSys_MassPoint);
  m.def("cppFEM_Merge_PointContact",      &PyMergeLinSys_Contact);
  m.def("cppFEM_Merge_ScalarPoission",    &PyMergeLinSys_Poission);
  m.def("cppFEM_Merge_ScalarDiffuse",     &PyMergeLinSys_Diffuse);
  m.def("cppFEM_Merge_SolidLinearStatic", &PyMergeLinSys_LinearSolidStatic);
  m.def("cppFEM_Merge_SolidLinearDynamic",&PyMergeLinSys_LinearSolidDynamic);
  m.def("cppFEM_Merge_FluidStorksStatic", &PyMergeLinSys_StorksStatic2D);
  m.def("cppFEM_Merge_FluidStorksDynamic",&PyMergeLinSys_StorksDynamic2D);
  m.def("cppFEM_Merge_FluidNavierStorks", &PyMergeLinSys_NavierStorks2D);
  m.def("cppFEM_Merge_ShellCloth",        &PyMergeLinSys_Cloth);
  m.def("cppFEM_Merge_ShellMitc3Static",  &PyMergeLinSys_ShellMitc3Static);

  m.def("pbd_proj_rigid2d",            &PyPBD_ConstProj_Rigid2D);
  m.def("pbd_proj_rigid3d",            &PyConstProj_Rigid3D);
  m.def("pbd_proj_cloth_stretch",      &PyPBD_ConstProj_ClothStretch);
  m.def("pbd_proj_cloth_bend",         &PyPBD_ConstProj_ClothBend);
  m.def("pbd_proj_seam",               &PyPBD_ConstProj_Seam);
  m.def("pbd_proj_contact",            &PyPBD_ConstProj_Contact);
  m.def("pbd_pointFixBC",              &PyPointFixBC);
}
