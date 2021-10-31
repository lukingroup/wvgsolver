from abc import ABC, abstractmethod

class Parser(ABC):
  def __init__(self, data, meta=None):
    self.data = data
    self.meta = meta

    self._parse()
  
  def _parse(self):
    pass

  @abstractmethod
  def show(self, **kwargs):
    pass

  @abstractmethod
  def save(self, fpath, **kwargs):
    pass

class EngineSpecificParser(Parser, ABC):
  def __init__(self, data, meta=None, engine=None):
    self.engine = engine
    super().__init__(data, meta)

  def set_engine(self, engine):
    self.engine = engine