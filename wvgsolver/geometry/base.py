from ..utils.base import EngineSpecific
from ..utils.misc import randstring
from abc import ABC, abstractmethod
from ..utils.constants import C_LIGHT
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
    super().__init__()
    self.order = order
    self.color = color
  
  def __repr__(self):
    return "Material(%d,%s):%s" % (self.order, self.color, super().__repr__())
