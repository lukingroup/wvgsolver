from abc import ABC
from .misc import hasmethod

class EngineSpecific(ABC):
  def callImplementation(self, name, session, *args, **kwargs):
    engine_name = session.engine.name
    fname = "_" + name + "_" + engine_name
    if not hasmethod(self, fname):
      raise NotImplementedError("Function '%s' is not implemented for engine '%s'" % \
              (name, engine_name))

    return getattr(self, fname)(session, *args, **kwargs)
