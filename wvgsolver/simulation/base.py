import uuid
from abc import ABC, abstractmethod
from ..utils.misc import hasmethod
from ..engine import getDefaultEngine

class SimulationObject(ABC):
  def __init__(self, engine=None):
    if engine is None:
      engine = getDefaultEngine()

    self.engine = engine
    self.results = None
    self.default_sim_type = ""

  @abstractmethod
  def _get_structures(self):
    return []

  def simulate(self, t=None, engine_options={}, **kwargs):
    if t is None:
      t = self.default_sim_type

    sim_type = "_simulate_" + t
    if not hasmethod(self, sim_type):
      raise NotImplementedError("simulation type '%s' not found for object '%s'" % \
              (t, type(self).__name__))

    session = self.engine.new_session()

    try:
      session.set_structures(self._get_structures())
      session.set_engine_options(engine_options)

      self.results = getattr(self, sim_type)(session, **kwargs)
    except Exception:
      session.close()
      raise

    session.close()

    return self.results

  def clear(self):
    self.results = None

