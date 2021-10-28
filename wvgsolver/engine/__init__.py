from .engines import DummyEngine, LumericalEngine

DEFAULT_ENGINE = DummyEngine()

# TODO: Check that this doesn't instantiate many copies of LumericalEngine
def getDefaultEngine():
  return DEFAULT_ENGINE

def setDefaultEngine(engine):
  global DEFAULT_ENGINE
  DEFAULT_ENGINE = engine

