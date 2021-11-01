from ..utils.base import EngineSpecific
from ..utils.misc import randstring
from abc import ABC, abstractmethod
from ..utils.constants import C_LIGHT
from matplotlib.colors import to_rgba
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
    self.callImplementation("add", session)

class Structure(Geometry, ABC):
  def __init__(self, pos, material):
    super().__init__()
    self.pos = pos
    self.material = material

  @abstractmethod
  def get_mesh(self, scale):
    """
    Gets this polygon's 3D mesh.

    Parameters
    ----------
    scale : number
      What size to treat as unity in the scale of the structure. For example, a value of 1e-6
      corresponds to a 1 micron scaling.

    Returns
    -------
    mesh : trimesh.Trimesh
    """
    pass

  def __repr__(self):
    return "Structure(%s):%s" % (self.pos, super().__repr__())

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
