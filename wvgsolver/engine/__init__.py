from .engines import LumericalEngine

DEFAULT_ENGINE = LumericalEngine()

def getDefaultEngine():
  return DEFAULT_ENGINE

def setDefaultEngine(engine):
  global DEFAULT_ENGINE
  DEFAULT_ENGINE = engine

