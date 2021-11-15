from .engines import DummyEngine, LumericalEngine
from .misc import LumericalMeshRegion

DEFAULT_ENGINE = DummyEngine()

# TODO: Check that this doesn't instantiate many copies of LumericalEngine
def getDefaultEngine():
  return DEFAULT_ENGINE

def setDefaultEngine(engine):
  global DEFAULT_ENGINE
  DEFAULT_ENGINE = engine

