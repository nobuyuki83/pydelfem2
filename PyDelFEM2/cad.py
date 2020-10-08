####################################################################
# Copyright (c) 2020 Nobuyuki Umetani                              #
#                                                                  #
# This source code is licensed under the MIT license found in the  #
# LICENSE file in the root directory of this source tree.          #
####################################################################

import numpy, os
from typing import Tuple, List

from .c_core import CppCad2D, CppMeshDynTri2D, CppMesher_Cad2D, CppVoxelGrid, CppMapper, AABB3
from .c_core import CAD_EDGE_GEOM_LINE, CAD_EDGE_GEOM_BEZIER_CUBIC, CAD_EDGE_GEOM_BEZIER_QUADRATIC
from .c_core import cppCad2D_ImportSVG, cppSVG_Polyline
from .c_core import cad_getPointsEdge
from .c_core import TRI, QUAD, HEX, TET, LINE
'''
from .c_core import \
  meshquad3d_voxelgrid, \
  meshquad3d_subdiv, \
  meshhex3d_voxelgrid, \
  meshhex3d_subdiv,\
  meshdyntri2d_initialize
from .c_core import meshtri3d_read_ply, meshtri3d_read_obj, meshtri3d_read_nastran, meshtri3d_write_obj
from .c_core import setXY_MeshDynTri2D
from .c_core import cad_getPointsEdge, cppJArray_MeshPsup, quality_meshTri2D
from .c_core import copyMeshDynTri2D
from .c_core import setTopology_ExtrudeTri2Tet
from .c_core import cppNormalVtx_Mesh, cppEdge_Mesh
'''

from .c_core import cppMvc
from .c_core import numpyXYTri_MeshDynTri2D

from .msh import MeshDynTri2D


####################

class Cad2D():
  def __init__(self):
    self.ccad = CppCad2D()

  def draw(self) -> None:
    self.ccad.draw()

  def mouse(self,btn,action,mods,src,dir,view_height) -> None:
    if btn == 0:
      if action == 1:
        self.ccad.pick(src[0],src[1],view_height)

  def motion(self,src0,src1,dir) -> None:
    self.ccad.drag_picked(src1[0],src1[1], src0[0],src0[1])

  def minmax_xyz(self):
    return self.ccad.minmax_xyz();

  ######

  def clear(self) -> None:
    """clear all the cad elements"""
    self.ccad.clear()

  def pick(self, x, y, view_height) -> None:
    self.ccad.pick(x,y,view_height)

  def add_polygon(self,list_xy) -> None:
    self.ccad.add_polygon(list_xy)
    self.ccad.check()

  def add_vtx_edge(self, iedge, pos:List[float]) -> None:
    self.ccad.add_vtx_edge(pos[0],pos[1],iedge)
    self.ccad.check()

  def add_vtx_face(self, iface, pos:List[float]) -> None:
    self.ccad.add_vtx_face(pos[0],pos[1],iface)
    self.ccad.check()

  def set_edge_type(self, iedge:int, type:int, param:List[float]):
    self.ccad.set_edge_type(iedge,type,param)

  def edge_type(self, iedge:int) -> int:
    return self.ccad.edge_type(iedge)

  def iedge_picked(self) -> int:
    return self.ccad.iedge_picked

  def ivtx_picked(self) -> int:
    return self.ccad.ivtx_picked

  def iface_picked(self) -> int:
    return self.ccad.iface_picked

  def clean_picked(self) -> None:
    self.ccad.ivtx_picked = -1
    self.ccad.iedge_picked = -1
    self.ccad.iface_picked = -1

  def points_edge(self, list_edge_index, np_xy, tolerance=0.01):
    return cad_getPointsEdge(self.ccad,list_edge_index, np_xy, tolerance=tolerance)

  def import_svg(self,path0:str,scale=(1.0,1.0)):
    if not os.path.isfile(path0):
      return false
    cppCad2D_ImportSVG(self.ccad,path0,scale[0],scale[1])
    self.ccad.check()

  def export_svg(self,path0:str,scale=1.0):
    list_xy = self.ccad.xy_vtxctrl_face(0)
    str0 = cppSVG_Polyline(list_xy,scale)
    with open(path0, mode='w') as f:
      f.write(str0)


########################################################################################

class CadMesh2D(Cad2D):
  def __init__(self,edge_length:float):
    super().__init__()
    self.ccad.is_draw_face = False
    self.edge_length = edge_length
    self.dmsh = MeshDynTri2D() # this object does not reallocate
    self.map_cad2msh = None # this object reallocate
    self.listW = list()
    self.is_sync_mesh = True
    self.mesher = Mesher_Cad2D(edge_length=edge_length)

  def draw(self):
    self.ccad.draw()
    self.dmsh.draw()

  def motion(self,src0,src1,dir):
    self.drag_picked(src1[0],src1[1], src0[0],src0[1])

  def minmax_xyz(self):
    return self.dmsh.minmax_xyz()

  #####

  def drag_picked(self, s1x,s1y, s0x,s0y):
    self.ccad.drag_picked(s1x,s1y, s0x,s0y)
    if not self.is_sync_mesh:
      return
    assert len(self.listW) == self.ccad.nface()
    for iface in range(self.ccad.nface()):
      list_xy_bound = self.ccad.xy_vtxctrl_face(iface)
      np_xy_bound = numpy.array(list_xy_bound).reshape([-1, 2])
      np_pos_face = numpy.dot(self.listW[iface][1],np_xy_bound)
      self.dmsh.np_pos[self.listW[iface][0]] = np_pos_face
      self.dmsh.syncXY_from_npPos()
      '''
      max_asp,min_area = quality_meshTri2D(self.dmsh.np_pos,self.dmsh.np_elm)
      if max_asp > 5.0 or min_area < 0.0:
        self.remesh()
      '''

  def remesh(self):
    self.mesher.meshing(self,self.dmsh)
    ####
    self.listW.clear()
    for iface in range(self.ccad.nface()):
      npIndPoint_face = self.mesher.points_on_faces([iface],self)
      npPosPoint_face = self.dmsh.np_pos[npIndPoint_face]
      np_xy_bound = numpy.array(self.ccad.xy_vtxctrl_face(iface)).reshape([-1, 2])
      W = cppMvc(npPosPoint_face, np_xy_bound)
      assert W.shape[0] == npPosPoint_face.shape[0]
      assert W.shape[1] == np_xy_bound.shape[0]
      self.listW.append( [npIndPoint_face,W] )
    assert len(self.listW) == self.ccad.nface()

  def add_vtx_edge(self, iedge:int, pos:List[float]):
    super().add_vtx_edge(iedge,[pos[0],pos[1]])
    self.remesh()

#  def add_polygon(self,list_xy):
#    self.ccad.add_polygon(list_xy)

#  def set_edge_type(self, iedge:int, type:int, param:List[float]):
#    super().set_edge_type(iedge,type,param)

#####################################################


class Mesher_Cad2D():
  def __init__(self,edge_length=0.01):
    self.cmshr = CppMesher_Cad2D()
    self.cmshr.edge_length = edge_length

  def points_on_faces(self,list_iface:List[int],cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_faces(list_iface,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def points_on_edges(self,list_iedge:List[int],cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_edges(list_iedge,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def points_on_one_edge(self,iedge:int,is_endpoints:bool, cad:Cad2D) -> numpy.ndarray:
    list_points = self.cmshr.points_on_one_edge(iedge,is_endpoints,cad.ccad)
    return numpy.array(list_points,dtype=numpy.int32)

  def meshing(self,cad:Cad2D,dmesh=None):
    cdmsh = CppMeshDynTri2D()
    self.cmshr.meshing(cdmsh,cad.ccad)
    np_pos, np_elm = numpyXYTri_MeshDynTri2D(cdmsh)
    if dmesh is None:
      dmesh = MeshDynTri2D()
    dmesh.cdmsh = cdmsh
    dmesh.np_pos = np_pos
    dmesh.np_elm = np_elm
    dmesh.elem_type = TRI
    return dmesh

