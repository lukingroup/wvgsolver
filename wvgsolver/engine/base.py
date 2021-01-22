from abc import ABC, abstractmethod
from ..utils.misc import randstring
import logging

class Engine(ABC):
  def __init__(self):
    self.name = "abstract"

  @abstractmethod
  def new_session(self):
    return None

class Session(ABC):
  def __init__(self, engine):
    self.name = randstring()
    self.engine = engine
    self.engine_options = {}
    self.analyzers = {}

    logging.debug("Started '%s' session with name '%s'" % (self.engine.name, self.name))

  def __del__(self):
    self.close()

  def set_engine_options(self, engine_options):
    self.engine_options = engine_options

  @abstractmethod
  def close(self):
    pass

  def set_structures(self, structs=[]):
    if isinstance(structs, list):
      return self._set_structures(structs)
    
    return self._set_structures([structs])
  
  @abstractmethod
  def _set_structures(self, structs=[]):
    pass

  @abstractmethod
  def set_sim_region(self, pos=None, size=None, boundaries={}):
    pass

  # Time in femtoseconds 
  @abstractmethod
  def set_sim_time(self, t):
    pass
  
  def set_sources(self, sources=[]):
    if isinstance(sources, list):
      return self._set_sources(sources)
    
    return self._set_sources([sources])
  
  @abstractmethod
  def _set_sources(self, sources=[]):
    pass

  @abstractmethod
  def _runsim(self, options={}):
    pass
  
  @abstractmethod
  def _prerun(self):
    pass

  def analyze(self, k, *args, **kwargs):
    if k not in self.analyzers:
      raise ValueError("Invalid analyze key '%s'" % k)

    return self.analyzers[k].analyze(self, *args, **kwargs)

  def run(self, analysis={}, analyze=True):
    self._prerun()

    if isinstance(analysis, dict):
      for a in analysis.values():
        a.setup(self)
    else:
      analysis.setup(self)

    self._runsim(self.engine_options)

    output = None
    if analyze:
      output = {}
      if isinstance(analysis, dict):
        for k, a in analysis.items():
          output[k] = a.analyze(self)
      else:
        output = analysis.analyze(self)
    else:
      self.analyzers = analysis

    return output

