from .base import Structure
import numpy as np
import shapely
import trimesh

class PolygonStructure(Structure):
  """Represents a 2D polygon extruded along the z axis into a 3D structure, then rotated and placed at a particular location.
  This allows for arbitrarily positioned extruded polygons, and is the most generic Structure as of right now
  """
  def __init__(self, pos, verts, height, material, rot_angles=(0, 0, 0)):
    """
    Parameters
    ----------
    pos : Vec3
      The position of the centroid of the structure, halfway up the extrusion 
    verts : array-like
      A n x 2 array of 2D points, corresponding to the x-y coordinates of the 2D polygonal face's
      vertices
    height : float
      The length of the extrusion of the polygonal face along the z-axis
    material : Material
      The material that the structure is made of
    rot_angles : tuple
      Euler angles in radians, given in extrinsic (static) z-y-x order, with which to rotate
      the 3D structure.  For example, a value of (pi/2, pi/2, 0) means rotate around the z axis
      by 90 degrees, then around the y axis by 90 degrees.
      When rotating, the origin in the x-y plane is taken to be the origin used when defining the vertices
      of the face, and the origin in the z dimension is halfway up the extrusion height of the structure.
    """
    super().__init__(pos, material)
    self.verts = verts
    self.rot_angles = rot_angles
    self.height = height

  def get_mesh(self, scale):
    transform = np.dot(
      trimesh.transformations.translation_matrix((self.pos / scale).tolist()),
      np.dot(
        trimesh.transformations.euler_matrix(*self.rot_angles, axes="szyx"),
        trimesh.transformations.translation_matrix([0, 0, -self.height / (2 * scale)])
      )
    )
    return trimesh.primitives.Extrusion(polygon=shapely.geometry.Polygon(np.array(self.verts) / scale), 
      height=(self.height / scale), transform=transform).to_mesh()

  def __repr__(self): 
    return "PolygonStructure(%d,%s,%.6e,%s):%s" % (len(self.verts), self.rot_angles, self.height, self.material, super().__repr__())

  def _add_lumerical(self, session):
    self.material.add(session)
    poly = session.fdtd.addpoly(name=self.name, x=self.pos.x, y=self.pos.y, z=self.pos.z,
      z_span=self.height, first_axis=4, second_axis=3, third_axis=2, material=self.material.name)
    poly.rotation_1 = self.rot_angles[0] * 180 / np.pi
    poly.rotation_2 = self.rot_angles[1] * 180 / np.pi
    poly.rotation_3 = self.rot_angles[2] * 180 / np.pi
    poly.vertices = np.array(self.verts)

class BoxStructure(PolygonStructure):
  """Represents a 3D box structure by building a Polygon structure from a 2D polygon defined by the 
  x-y dimensions of the box and an extrusion height defined by the z dimension of the box"""
  def __init__(self, pos, size, material, rot_angles=(0, 0, 0)):
    """
    Parameters
    ----------
    pos : Vec3
      The position of the center of the box
    size : Vec3
      The length of each side of the box
    material : Material
      The material that the structure is made of
    rot_angles : tuple
      See documentation of PolygonStructure
    """
    verts = [
      [-size.x/2, -size.y/2],
      [-size.x/2, size.y/2],
      [size.x/2, size.y/2],
      [size.x/2, -size.y/2]
    ]
    super().__init__(pos, verts, size.z, material, rot_angles)
    self.size = size

  def __repr__(self):
    return "BoxStructure(%s):%s" % (self.size, super().__repr__())

class CylinderStructure(PolygonStructure):
  """Represents an elliptical cylinder structure. Internally this is still a PolygonStructure,
  with a high vertex count polygon used to approximate the ellipse. When not rotated, this 
  cylinder has its axis oriented along the z axis"""
  def __init__(self, pos, height, radius, material, radius2=None, rot_angles=(0, 0, 0), ncirclepoints=100):
    """
    Parameters
    ----------
    pos : Vec3
      The position of the center of the cylinder
    height : float
      The height of the cylinder
    radius : float
      The radius of the cylinder along the x axis (before rotating)
    material : Material
      The material that the structure is made of
    radius2 : float or None
      If provided, the radius of the cylinder along the y axis (before rotating). If not provided,
      this is taken to be equal to radius
    rot_angles : tuple
      See documentation of PolygonStructure
    ncirclepoints : int
      The number of vertices to use in the polygon that approximates the elliptical face of the cylinder
    """
    self.ncirclepoints = ncirclepoints
    self.radius = radius
    self.radius2 = radius2 if radius2 is not None else radius
    verts = self._get_verts()
    
    super().__init__(pos, verts, height, material, rot_angles)

  def _get_verts(self):
    verts = []
    for i in range(self.ncirclepoints):
      t = (i / self.ncirclepoints) * 2 * np.pi
      verts.append([ self.radius * np.cos(t), self.radius2 * np.sin(t) ])
    return verts

  def __repr__(self):
    return "CylinderStructure(%.6e,%.6e, %d):%s" % (self.radius, self.radius2, self.ncirclepoints, super().__repr__())

  def _add_lumerical(self, session):
    self.material.add(session)
    circle = session.fdtd.addcircle(name=self.name, x=self.pos.x, y=self.pos.y, z=self.pos.z,
      z_span=self.height, radius=self.radius, material=self.material.name, make_ellipsoid=True,
      first_axis=4, second_axis=3, third_axis=2)
    circle.radius_2 = self.radius2
    circle.rotation_1 = self.rot_angles[0] * 180 / np.pi
    circle.rotation_2 = self.rot_angles[1] * 180 / np.pi
    circle.rotation_3 = self.rot_angles[2] * 180 / np.pi
