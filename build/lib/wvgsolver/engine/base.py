from abc import ABC, abstractmethod
from ..utils.misc import randstring
import logging
import copy

class Engine(ABC):
  def __init__(self):
    self.name = "abstract"

  def __repr__(self):
    return "Engine(%s)" % self.name

  @abstractmethod
  def new_session(self):
    return None

class Session(ABC):
  def __init__(self, engine):
    self.name = randstring()
    self.engine = engine
    self.analyzers = {}

    logging.info("Started '%s' session with name '%s'" % (self.engine.name, self.name))

  def __del__(self):
    self.close()

  @abstractmethod
  def close(self):
    pass

  def set_structures(self, structs=[]):
    if isinstance(structs, list):
      return self._set_structures(structs)
    
    return self._set_structures([structs])

  def set_mesh_regions(self, regions=[]):
    if isinstance(regions, list):
      return self._set_mesh_regions(regions)
    
    return self._set_mesh_regions([regions])
  
  @abstractmethod
  def _set_mesh_regions(self, structs=[]):
    pass

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

  def get_postrunres(self):
    return None

  def run(self, analysis={}, analyze=True):
    self._prerun()

    # Clean up the last run's analyzers
    if isinstance(self.analyzers, dict):
      for a in self.analyzers.values():
        a.cleanup(self)
    elif self.analyzers:
      self.analyzers.cleanup(self)
    self.analyzers = analysis

    if isinstance(analysis, dict):
      for a in analysis.values():
        a.setup(self)
    else:
      analysis.setup(self)

    self._runsim()

    output = None
    if analyze:
      output = {}
      if isinstance(self.analyzers, dict):
        for k, a in self.analyzers.items():
          output[k] = a.analyze(self)
      else:
        output = self.analyzers.analyze(self)

    return output

class SessionResult(ABC):
  def get_clean(self):
    engine = self.engine
    del self.engine
    c = copy.deepcopy(self)
    self.engine = engine
    return c
    
