from ..utils.base import EngineSpecific

class Analysis(EngineSpecific):
  def setup(self, session):
    self.callImplementation("setup", session, session)
  
  def analyze(self, session, *args, **kwargs):
    return self.callImplementation("analyze", session, session, *args, **kwargs)

