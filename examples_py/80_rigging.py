import PyDelFEM2 as dfm2
import PyDelFEM2.gl.glfw

def test_gltf():
  gltf = dfm2.CppGLTF()
  #gltf.read("../test_inputs/RiggedFigure.glb")
  gltf.read("../test_inputs/CesiumMan.glb")

  np_pos0,np_elm,np_rigw,np_rigj = dfm2.cppGetMeshInfoGltf(gltf,0,0)
  np_pos = np_pos0.copy()
  msh = dfm2.Mesh(np_pos,np_elm,dfm2.TRI)
  dfm2.gl.glfw.winDraw3d([msh],winsize=(400,300))

  skeleton = dfm2.cppGetSkeleton_Gltf(gltf,0)
#  skeleton.set_translation( 0, [0.0, 0.0, +0.2] )
  skeleton.set_rotation_bryant(0, [-3.1415*0.5, 0.0, -3.1415*0.5] )
#  skeleton.set_rotation_bryant(2, [-3.1415*0.0, 0.0, 0.0] )
  skeleton.set_rotation_bryant(5, [0.5, 0.0, 0.0] )
  dfm2.cppUpdateRigSkin(np_pos,
                        np_pos0, skeleton, np_rigw, np_rigj)
  msh = dfm2.Mesh(np_pos,np_elm,dfm2.TRI)
  axis = dfm2.gl.AxisXYZ(1)
  dfm2.gl.glfw.winDraw3d([msh,axis],winsize=(400,300))


def test_smpl0():
  import numpy as np
  smpl = np.load("../test_inputs/smpl_model_f.npz")
  npV = smpl["vertices_template"]
  npT = smpl["face_indices"]
  npT[:,:] -= 1 # index should be minus 1
  msh = dfm2.Mesh(npV,npT,dfm2.TRI)
  axis = dfm2.gl.AxisXYZ(1)
  dfm2.gl.glfw.winDraw3d([msh,axis],winsize=(400,300))

  npIdBoneParent = smpl["kinematic_tree"]
  npJ = np.matmul(smpl["joint_regressor"],npV) # joint position
  npSparseW,npSparseI = dfm2.cppSparsifyMatrixRow(smpl["weights"])

  skeleton = dfm2.cppGetSkeleton_BoneTreePos(npIdBoneParent[0],npJ)
  #skeleton.set_translation( 0, [0.0, 0.0, +0.2] )
  skeleton.set_rotation_bryant(14, [0.0, 0.0, +3.1415*0.3] )
  npV1 = npV.copy()
  dfm2.cppUpdateRigSkin(npV1,
                        npV, skeleton, npSparseW, npSparseI)
  msh = dfm2.Mesh(npV1,npT,dfm2.TRI)
  axis = dfm2.gl.AxisXYZ(1)
  dfm2.gl.glfw.winDraw3d([msh,axis],winsize=(400,300))


def test_smpl1():
  import numpy as np
  smpl = np.load("../test_inputs/smpl_model_f.npz")
  npV = smpl["vertices_template"]
  npT = smpl["face_indices"]
  npT[:,:] -= 1 # index should be minus 1
  npIdBoneParent = smpl["kinematic_tree"]
  npJntRgrs = smpl["joint_regressor"]
  print(npJntRgrs.shape)
  npSparseRigW, npSparseRigI = dfm2.cppSparsifyMatrixRow(smpl["weights"])
  npSparseJntW, npSparseJntI = dfm2.cppSparsifyMatrixRow(smpl["joint_regressor"])
  print("npSparseRig",npSparseRigW.shape, npSparseRigI.shape)
  print("npSparseJnt",npSparseJntW.shape, npSparseJntI.shape)

  npBlendShape = smpl["shape_blend_shapes"]
  npBlendPose = smpl["pose_blend_shapes"]





if __name__ == "__main__":
  test_gltf()
  test_smpl0()
  test_smpl1()