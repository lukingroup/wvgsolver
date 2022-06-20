from .misc import DummyEngine
from .lumerical import LumericalEngine
from .meep import EffIndex1DEngine

DEFAULT_ENGINE = DummyEngine()

# TODO: Check that this doesn't instantiate many copies of LumericalEngine
def getDefaultEngine():
  return DEFAULT_ENGINE

def setDefaultEngine(engine):
  global DEFAULT_ENGINE
  DEFAULT_ENGINE = engine

