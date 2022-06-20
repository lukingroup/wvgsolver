from .base import Engine

class DummyEngine(Engine):
  """An engine that does nothing"""
  def __init__(self):
    self.name = "dummy"

  def new_session(self):
    raise NotImplementedError("Dummy engine doesn't support creating sessions. Please use LumericalEngine or other")
