from ..utils.base import EngineSpecific
from ..utils.misc import randstring
from abc import ABC, abstractmethod
from ..utils.constants import C_LIGHT
from matplotlib.colors import to_rgba
import trimesh
import numpy as np
import copy

class Geometry(EngineSpecific):
  def __init__(self):
    self.name = randstring()

  def copy(self):
    g = copy.deepcopy(self)
    g.name = randstring()
    g.material.name = randstring()
    return g

  def __repr__(self):
    return "Geometry(%s)" % self.name

  def add(self, session):
    return self.callImplementation("add", session)

class Structure(Geometry, ABC):
  def __init__(self, pos, material, rot_angles=(0, 0, 0)):
    """
    Parameters
    ----------
    pos : Vec3
      The position of the centroid of the structure
    material : Material
      The material that the structure is made of
    rot_angles : tuple
      Euler angles in radians, given in extrinsic (static) z-y-x order, with which to rotate
      the 3D structure.  For example, a value of (pi/2, pi/2, 0) means rotate around the z axis
      by 90 degrees, then around the y axis by 90 degrees.
    """
    self.pos = pos
    self.material = material
    self.rot_angles = rot_angles
    
    super().__init__()

  def get_mesh(self, scale):
    """
    Gets this structure's 3D mesh.

    Parameters
    ----------
    scale : number
      What size to treat as unity in the scale of the structure. For example, a value of 1e-6
      corresponds to a 1 micron scaling.

    Returns
    -------
    mesh : trimesh.Trimesh
    """
    transform = np.dot(
      trimesh.transformations.translation_matrix((self.pos / scale).tolist()),
      trimesh.transformations.euler_matrix(*self.rot_angles, axes="szyx"),
    )

    mesh = self._get_origin_mesh(scale)
    mesh.apply_transform(transform)
    mesh.visual = trimesh.visual.ColorVisuals(
      face_colors=np.tile(
        np.minimum(255, np.floor(256*self.material.get_color_rgba())).astype(np.uint8),
        (len(mesh.faces), 1)
      )
    )
    return mesh

  @abstractmethod
  def _get_origin_mesh(self, scale):
    """
    Returns the mesh of this object, centered at the origin, and not rotated

    Returns
    -------
    mesh : trimesh.Trimesh
    """
    pass

  @abstractmethod
  def contains(self, points, scale=1e-6):
    pass

  def __repr__(self):
    return "Structure(%s):%s" % (self.pos, super().__repr__())
  
  def __hash__(self):
    return hash(str(self.pos) + str(hash(self.material)) + str(self.rot_angles) + \
      str(trimesh.comparison.identifier_simple(self._get_origin_mesh(1e-6))))

  def __eq__(self, other):
    if not isinstance(other, Structure):
      return False

    if self.pos != other.pos or self.material != other.material or self.rot_angles != other.rot_angles:
      return False

    comp = trimesh.comparison.identifier_simple

    scale = 1e-6
    mesh = comp(self._get_origin_mesh(1e-6))
    other_mesh = comp(other._get_origin_mesh(1e-6))

    return np.array_equal(mesh, other_mesh)

class Source(Geometry, ABC):
  def __init__(self, pos, frange, f, pulse_length, pulse_offset):
    super().__init__()
    self.pos = pos
    self.frange = frange
    self.f = f
    self.pulse_length = pulse_length
    self.pulse_offset = pulse_offset

    if frange is not None:
      self.lambda_start = C_LIGHT / frange[1]
      self.lambda_end = C_LIGHT / frange[0]
  
  def _config_freq_lumerical(self, source):
    source.name = self.name
    source.x = self.pos.x
    source.y = self.pos.y
    source.z = self.pos.z

    if self.frange is not None:
      source.wavelength_start = C_LIGHT / self.frange[1]
      source.wavelength_stop = C_LIGHT / self.frange[0]
    else:
      source.set_time_domain = True
      source.pulse_type = 1
      source.frequency = self.f
      if self.pulse_length is not None:
        source.pulselength = self.pulse_length
      if self.pulse_offset is not None:
        source.offset = self.pulse_offset

class Material(Geometry, ABC):
  def __init__(self, order, color):
    """
    Parameters
    ----------
    order : number
      The override priority of the material. Must be an integer, greater than or equal to 0.
      A lower order corresponds to a higher priority, so if there exists a region where two 
      structures overlap the structure with the lower order material will remain. 
    color : Matplotlib color
      The color of material when visualized. Can be a string "blue", "red", etc. or any valid
      matplotlib color (see the matplotlib docs for more info).
    """ 

    self.order = order
    self.color = color
    super().__init__()
 
  def get_color_rgba(self):
    """Returns the color as an RGBA array, with each value from 0 to 1"""
    return np.array(to_rgba(self.color))
  
  def __repr__(self):
    return "Material(%d,%s):%s" % (self.order, self.color, super().__repr__())

  def __hash__(self):
    return hash(self.order)

  def __eq__(self, other):
    return self.order == other.order
