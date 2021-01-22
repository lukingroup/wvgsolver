from .base import Engine
from .sessions import LumericalSession
import importlib.util
import logging
import sys
import os
import tempfile

class LumericalEngine(Engine):
  def __init__(self, lumerical_path=None, working_path=None):
    self.name = "lumerical"
    self.lumapi = None
    self.lumapi_path = lumerical_path
    if working_path is None:
      self.working_path = os.path.join(tempfile.gettempdir(), "wvgsolver")
    else:
      self.working_path = working_path

    if not os.path.isdir(self.working_path):
      os.mkdir(self.working_path)

  def checkapi(self):
    if self.lumapi_path is None:
      base_path = ""
      if sys.platform == "linux":
        base_path = "/opt/lumerical/"
      elif sys.platform == "win32" or sys.platform == "cygwin":
        base_path = "C:\\Program Files\\Lumerical"
      elif sys.platform == "darwin":
        self.lumapi_path = "/Applications/Lumerical/FDTD Solutions/FDTD Solutions.app/Contents/API/Python/lumapi.py"

      if base_path != "":
        if not os.path.isdir(base_path):
          raise RuntimeError("Could not find lumerical installation in default path")

        versions = sorted([
          d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and not d.startswith(".")
        ])

        if len(versions) < 1:
          raise RuntimeError("Could not find lumerical installation in default path")

        logging.debug("Using lumerical version '%s'" % versions[0])

        self.lumapi_path = os.path.join(base_path, versions[0], "api", "python", "lumapi.py")
      
      if not os.path.isfile(self.lumapi_path):
        raise RuntimeError("Invalid lumerical api path '%s', file not found" % self.lumapi_path)

      if self.lumapi is None:
        spec = importlib.util.spec_from_file_location("lumapi", self.lumapi_path)
        self.lumapi = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.lumapi)

  def new_session(self):
    self.checkapi()
    return LumericalSession(self)

