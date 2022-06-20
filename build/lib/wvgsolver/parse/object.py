from .base import Parser
from ..utils.linalg import Vec3
from ..utils.misc import randstring
from ..utils.constants import AXIS_X, AXIS_Y, AXIS_Z
from ..geometry.materials import DielectricMaterial
from ..geometry.structures import PolygonStructure
import gdspy
import trimesh
import numpy as np
import shapely
import logging
import phidl
from abc import ABC, abstractmethod
phidl.set_quickplot_options(blocking=True)

class ObjectParser(Parser, ABC):
  """Parses a SimlationObject for some kind of visualization task"""
  def __init__(self, data, scale=1e-6, origin=Vec3(0), layer=0):
    """
    Parameters
    ----------
    data : SimulationObject
      The simulation object being parsed
    scale : float
      The scale of the resulting GDS. A value of 1e-6 means an object of length 1 in the GDS 
      corresponds to geometry of length 1 micron
    origin : Vec3
      The origin of the slicing plane
    layer : int
      Layer of the resulting geometry
    """
    self.scale = scale
    self.origin = origin.tolist()
    self.layer = layer

    super().__init__(data)

class ObjectGDSParser(ObjectParser, ABC):
  """Base class for parsers that take a SimulationObject, slice it along a plane, and return 
  the resulting 2D geometry as a Phidl GDS device"""
  def __init__(self, data, scale=1e-6, layer=0, origin=Vec3(0), axis=AXIS_Z, name="main", invert=False):
    """
    Parameters
    ----------
    data : SimulationObject
      See documentation of ObjectParser
    scale : float
      See documentation of ObjectParser
    layer : int
      See documentation of ObjectParser
    origin : Vec3
      See documentation of ObjectParser
    axis : {AXIS_X, AXIS_Y, AXIS_Z}
      The axis of the normal of the slicing plane
    name : str
      The name of the Phidl device
    invert : bool
      Whether or not to invert the geometry at the end
    """
    self.invert = invert
    self.axis = axis
    self.device = phidl.Device(name=name)
    self.polygonset = None
    
    super().__init__(data, scale, origin, layer)

  def show(self):
    """Show the geometry in a plot"""
    phidl.quickplot(self.device)
  
  def save(self, fpath):
    """Save the geometry to a GDS file

    Parameters
    ----------
    fpath : str
      The file path to save to
    """
    self.device.write_gds(fpath)

class DielectricExtrusionFaceGDSParser(ObjectGDSParser):
  """
  Builds a GDS file containing all the extrusion faces of the structures in the given SimulationObject, on the same plane.
  """
  def _parse(self):
    structs = self.data.get_structures()

    polys = []
    origin = np.array(self.origin) / self.scale
    for s in structs:
      if not isinstance(s, PolygonStructure):
        logging.warn("Skipping object %s because it's not a polygon" % s)
        continue
      if not isinstance(s.material, DielectricMaterial):
        logging.warn("Skipping object %s because its material is not a dielectric" % s)
        continue
      polys.append({
        "poly": gdspy.Polygon(((np.array(s.verts) + np.array([s.pos.x, s.pos.y])) / self.scale) - origin[:2]),
        "struct": s
      })

    polys.sort(reverse=True, key=lambda p: p["struct"].material.order)

    if polys:
      totalpoly = gdspy.PolygonSet([])
      for pidx in range(len(polys)):
        op = "or" if polys[pidx]["struct"].material.nindex > 1  else "not"
        if self.invert:
          op = "not" if op == "or" else "or"
        totalpoly = gdspy.boolean(totalpoly, polys[pidx]["poly"], op, layer=self.layer)
      
      if totalpoly:
        self.polygonset = totalpoly
        self.device.add_polygon(totalpoly)

class DielectricConvexSliceGDSParser(ObjectGDSParser):
  """This class slices a Simulation Object as described in ObjectGDSParser, and returns the geometry resulting from 
  the intersection of all structures with an index of refraction greater than 1 with the plane. If invert=True, 
  the parser instead returns the geometry from the intersection of all air holes or structures with index of refraction 1 
  with the plane.

  NOTE: This parser assumes each individual structure in the given SimulationObject's geometry is convex under the necessary 
  2D slice. 
  """
  def _parse(self):
    axis_map = {}
    axis_map[AXIS_X] = 0
    axis_map[AXIS_Y] = 1
    axis_map[AXIS_Z] = 2

    structs = self.data.get_structures()

    polys = []
    origin = np.array(self.origin) / self.scale
    for s in structs:
      if not isinstance(s.material, DielectricMaterial):
        logging.warn("Skipping object %s because its material is not a dielectric" % s)
        continue
      idxs = [ i for i in range(3) if axis_map[self.axis] != i]
      intersection = (
        trimesh.intersections.mesh_plane(
            s.get_mesh(self.scale), 
            [int(self.axis == AXIS_X), int(self.axis == AXIS_Y), int(self.axis == AXIS_Z)], 
            origin
          )[:,:,idxs] - origin[idxs]
      ).tolist()
      hull = shapely.geometry.MultiLineString(intersection).convex_hull
      if not isinstance(hull, shapely.geometry.Polygon):
        logging.warn("Skipping object %s because its convex hull is not a polygon" % s)
        continue
      polys.append({
        "poly": gdspy.Polygon(list(hull.exterior.coords)[0:-1]),
        "struct": s
      })

    polys.sort(reverse=True, key=lambda p: p["struct"].material.order)

    if polys:
      totalpoly = gdspy.PolygonSet([])
      for pidx in range(len(polys)):
        op = "or" if polys[pidx]["struct"].material.nindex > 1  else "not"
        if self.invert:
          op = "not" if op == "or" else "or"
        totalpoly = gdspy.boolean(totalpoly, polys[pidx]["poly"], op, layer=self.layer)
      
      if totalpoly:
        self.polygonset = totalpoly
        self.device.add_polygon(totalpoly)

class ObjectParser3D(ObjectParser):
  """
  Takes a SimulationObject and renders it in 3D, either to a window or to an image.
  To avoid Z-fighting, this scales structures with lower order a tiny bit larger than
  structures with a higher order.
  """
  def __init__(self, data, scale=1e-6, origin=Vec3(0), order_scale=0.01):
    """
    Parameters
    ----------
    data : SimulationObject
      See documentation of ObjectParser
    scale : float
      See documentation of ObjectParser
    origin : Vec3
      See documentation of ObjectParser
    order_scale : float
      How much the order affects the scale (to prevent z-fighting)
    """

    self.scene = trimesh.scene.Scene()
    self.order_scale = order_scale
    super().__init__(data, scale, origin)

  def _parse(self):
    structs = self.data.get_structures()

    min_order = None
    for s in structs:
      if min_order is None or s.material.order < min_order:
        min_order = s.material.order

    for s in structs:
      mesh = s.get_mesh(self.scale * (1 + self.order_scale * (s.material.order - min_order)))
      self.scene.add_geometry(mesh)

  def show(self):
    """Show the object in a 3D visualizer"""
    self.scene.show()

  def save(self, fpath):
    # TODO: implement
    pass
