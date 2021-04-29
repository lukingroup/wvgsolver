from ..utils.base import EngineSpecific

class Analysis(EngineSpecific):
  def setup(self, session):
    self.callImplementation("setup", session)

  def analyze(self, session, *args, **kwargs):
    return self.callImplementation("analyze", session, *args, **kwargs)

  def cleanup(self, session):
    return self.callImplementation("cleanup", session)
