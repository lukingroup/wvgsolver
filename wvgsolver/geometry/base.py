from ..utils.base import EngineSpecific
from ..utils.misc import randstring
from abc import ABC
import copy

class Geometry(EngineSpecific):
  def __init__(self):
    self.name = randstring()

  def copy(self):
    g = copy.deepcopy(self)
    g.name = randstring()
    return g

  def add(self, session):
    self.callImplementation("add", session, session)

class Structure(Geometry, ABC):
  def __init__(self, pos):
    super().__init__()
    self.pos = pos

class Source(Geometry, ABC):
  def __init__(self, pos):
    super().__init__()
    self.pos = pos

class Shape(ABC):
  pass

class Material(Geometry, ABC):
  def __init__(self):
    super().__init__()
