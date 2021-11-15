from .base import Engine
from .sessions import LumericalSession
import importlib.util
import logging
import sys
import os
import tempfile
import subprocess
import time

class DummyEngine(Engine):
  """An engine that does nothing"""
  def __init__(self):
    self.name = "dummy"

  def new_session(self):
    raise NotImplementedError("Dummy engine doesn't support creating sessions. Please use LumericalEngine or other")

class LumericalEngine(Engine):
  """Engine for running simulations in Lumerical"""
  def __init__(self, lumerical_path=None, working_path=None, hide=True, save_fsp=True, mesh_accuracy=5, pml_layers=8):
    """
    Parameters
    ----------
    lumerical_path : str or None
      If provided, the file path to the root of your Lumerical installation. If not provided,
      we attempt to infer this from your operating system
    working_path : str or None
      If provided, the path where all FSP files in Analysis mode from simulations run 
      with this engine will be stored. If not provided, a temporary directory is created 
      for this, which is later deleted upon deletion of the Engine instance.
    hide : bool
      If True, don't display the Lumerical GUI while running simulations
    save_fsp : bool
      If True, save FSP files (in Layout mode) into memory so that they are accessible from the 
      list of simulation results in any SimulationObject that runs simulations using this engine. 
    mesh_accuracy : int
      The autmoeshing accuracy, from 1 to 8. Higher means a finer mesh.
    pml_layers : int
      The number of PML layers to use in PML boundary conditions
    """
    self.name = "lumerical"
    self.lumapi = None
    self.lumerical_path = lumerical_path
    self.hide = hide
    self.mesh_accuracy = mesh_accuracy
    self.pml_layers = pml_layers
    self.save_fsp = save_fsp
    self.version = "user specified" if lumerical_path is not None else "unknown"

    if working_path is None:
      self.temp_dir = tempfile.TemporaryDirectory()
      self.working_path = self.temp_dir.name
    else:
      self.temp_dir = None
      self.working_path = working_path
      if not os.path.isdir(self.working_path):
        os.mkdir(self.working_path)
    
    self.working_path = os.path.abspath(self.working_path)

  def __repr__(self):
    return "LumericalEngine(%s, %s, %s, %s, %d, %d, %s):%s" % (self.lumerical_path, self.working_path, self.hide, self.save_fsp,
      self.mesh_accuracy, self.pml_layers, self.version, super().__repr__())

  def open_fsp(self, fpath):
    self.checkapi()

    fdtdpath = os.path.join(self.lumerical_path, "bin", "fdtd-solutions")
    subprocess.run([fdtdpath, fpath], check=True)

  def checkapi(self):
    if self.lumapi is not None:
      return

    if self.lumerical_path is None:
      base_path = ""
      if sys.platform == "linux":
        base_path = "/opt/lumerical/"
      elif sys.platform == "win32" or sys.platform == "cygwin":
        base_path = "C:\\Program Files\\Lumerical"
      elif sys.platform == "darwin":
        self.lumerical_path = "/Applications/Lumerical/FDTD Solutions/FDTD Solutions.app/Contents"
        lumapi_path = os.path.join(self.lumerical_path, "API", "Python", "lumapi.py")

      if base_path != "":
        if not os.path.isdir(base_path):
          raise RuntimeError("Could not find lumerical installation in default path")

        versions = sorted([
          d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and not d.startswith(".")
        ])

        if len(versions) < 1:
          raise RuntimeError("Could not find lumerical installation in default path")
        version = versions[0]

        self.lumerical_path = os.path.join(base_path, version)
        logging.debug("Using lumerical version '%s'" % version)

        lumapi_path = os.path.join(self.lumerical_path, "api", "python", "lumapi.py")
        self.version = version
    else:
      lumapi_path = os.path.join(self.lumerical_path, "api", "python", "lumapi.py")

    if not os.path.isfile(lumapi_path):
      raise RuntimeError("Invalid lumerical api path '%s', file not found" % lumapi_path)

    spec = importlib.util.spec_from_file_location("lumapi", lumapi_path)
    self.lumapi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(self.lumapi)

  def new_session(self):
    self.checkapi()
    return LumericalSession(self)

