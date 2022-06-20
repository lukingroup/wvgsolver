from .base import Engine, Session
import importlib.util
import logging
import sys
import os
import tempfile
import subprocess
import time
import threading
import glob
from ..parse.engine import LumericalFSPFile

# TODO: Check why this is breaking, and make it delete 
# log files to be more robust
def lumericalLogFile(dir_path, name, stop):
  fps = {}
  while not stop.is_set():
    for f in glob.iglob(os.path.join(dir_path, name + "*.log")):
      try:
        base_name = os.path.basename(f)
        log_name = base_name[len(name):base_name.index(".log")]
        if log_name[0] == "_":
          log_name = log_name[1:]

        if f not in fps:
          fps[f] = open(f, "r")

        while True:
          line = fps[f].readline()
          if line:
            logging.info("[" + log_name + "] " + line.rstrip("\n\r"))
          else:
            break
      except Exception:
        logging.exception("Error reading log file '%s'" % f)

    time.sleep(1)

  for f, fp in fps.items():
    fp.close()
    if os.path.isfile(f):
      os.remove(f)

  return 0

class LumericalSession(Session):
  def __init__(self, engine):
    super().__init__(engine)
    self.fdtd = None
    self.sim_region = None

    self.structures_group = "structures"
    self.sources_group = "sources"
    self.mesh_regions_group = "mesh regions"
    self.boundary_keys_map = {
      "xmin": "x min bc",
      "xmax": "x max bc",
      "ymin": "y min bc",
      "ymax": "y max bc",
      "zmin": "z min bc",
      "zmax": "z max bc"
    }
    self.boundary_values_map = {
      "pml": 1,
      "metal": 2,
      "periodic": 3,
      "symmetric": 4,
      "antisymmetric": 5,
      "bloch": 6,
      "pmc": 7
    }

    if self.engine.temp_dir:
      self.temp_dir = tempfile.TemporaryDirectory(dir=self.engine.working_path)
      self.working_path = self.temp_dir.name
    else:
      self.temp_dir = None
      self.working_path = os.path.join(self.engine.working_path, self.name)
      if not os.path.isdir(self.working_path):
        os.mkdir(self.working_path)

    self.working_path = os.path.abspath(self.working_path)
    self.save_path = os.path.join(self.working_path, self.name + ".fsp")
    self.fsp_data = None

    self.fdtd = self.engine.lumapi.FDTD(hide=self.engine.hide, serverArgs={ "timeoutdasda": "0" })
    self.fdtd.addstructuregroup(name=self.structures_group)
    self.fdtd.addgroup(name=self.sources_group)
    self.fdtd.addgroup(name=self.mesh_regions_group)
    self.sim_region = self.fdtd.addfdtd(mesh_accuracy=self.engine.mesh_accuracy, use_early_shutoff=False)
    self.sim_region.pml_profile = 4
    self.sim_region.pml_min_layers = self.engine.pml_layers
    self.sim_region.pml_max_layers = self.engine.pml_layers

  def _set_mesh_regions(self, regions=[]):
    self.fdtd.switchtolayout()
    self.fdtd.groupscope(self.mesh_regions_group)
    self.fdtd.deleteall()

    for r in regions:
      r.add(self)
    
    self.fdtd.groupscope("::model")

  def close(self):
    if self.fdtd is not None:
      self.fdtd.close()

    self.fdtd = None
    self.sim_region = None
  
  def _set_structures(self, structs=[]):
    self.fdtd.switchtolayout()
    self.fdtd.groupscope(self.structures_group)
    self.fdtd.deleteall()

    for s in structs:
      s.add(self)
    
    self.fdtd.groupscope("::model")
  
  def _set_sources(self, sources=[]):
    self.fdtd.switchtolayout()
    self.fdtd.groupscope(self.sources_group)
    self.fdtd.deleteall()

    for s in sources:
      s.add(self)

    self.fdtd.groupscope("::model")

  def set_sim_time(self, t):
    self.fdtd.switchtolayout()
    self.sim_region.simulation_time = t

  def set_sim_region(self, pos=None, size=None, boundaries={}, dim2=False):
    self.fdtd.switchtolayout()
    if pos is not None:
      self.sim_region.x = pos.x
      self.sim_region.y = pos.y
      self.sim_region.z = pos.z

    if size is not None:
      self.sim_region.x_span = size.x
      self.sim_region.y_span = size.y
      self.sim_region.z_span = size.z

    for key in boundaries.keys():
      if len(key) == 1:
        val = boundaries[key] if not isinstance(boundaries[key], dict) else boundaries[key]["type"]
        self.sim_region[key + " min bc"] = self.boundary_values_map[val]
        if val == "bloch" and "k" in boundaries[key]:
          self.sim_region["set based on source angle"] = False
          self.sim_region["bloch units"] = 1
          self.sim_region["k" + key] = boundaries[key]["k"]
        elif val != "periodic":
          self.sim_region[key + " max bc"] = self.boundary_values_map[val]
      else:
        mapped_key = self.boundary_keys_map[key]
        min_key = mapped_key[0] + " min bc"
        if self.sim_region[min_key] in [3, 6]:
          self.sim_region[min_key] = 1
        self.sim_region[mapped_key] = self.boundary_values_map[boundaries[key]]
 
    self.sim_region.dimension = 1 if dim2 else 2

  def _prerun(self):
    self.fdtd.switchtolayout()
    for f in glob.iglob(os.path.join(self.working_path, self.name + "*.log")):
      try:
        os.remove(f)
      except:
        pass
  
  def _load_fsp(self):
    if self.engine.save_fsp:
      self.fsp_data = open(self.save_path, "rb").read()
    else:
      self.fsp_data = None

  def _runsim(self):
    self.fdtd.switchtolayout()
    self.fdtd.save(self.save_path)
    self._load_fsp()

#    stop_logging = threading.Event()
#    threading.Thread(target=lumericalLogFile, args=(self.working_path, self.name, stop_logging)).start()

    try:
      self.fdtd.run()
    except Exception:
#      stop_logging.set()
      raise

#    stop_logging.set()
  
  def get_postrunres(self):
    if self.fsp_data or not self.temp_dir:
      return LumericalFSPFile(self.fsp_data, self.save_path if not self.temp_dir else None, self.engine)
    return None

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

class LumericalMeshRegion:
  """Defines a mesh override region in Lumerical FDTD"""
  def __init__(self, bbox, dx, dy=None, dz=None):
    """
    Parameters
    ----------
    bbox : BBox
      The bbox of the region
    dx : float
      The maximum mesh step, in meters. To specificy the same step in all directions, 
      only provide this value. If dy and dz are specified, this is only taken to be
      the step size in the x direction.
    dy : float or None
      The maximum mesh step, in meters, in the y direction.
    dz : float or None
      The maximum mesh step, in meters, in the z direction.
    """
    self.bbox = bbox
    self.dx = dx
    self.dy = dx
    self.dz = dx
    if dy is not None:
      self.dy = dy
    if dz is not None:
      self.dz = dz 

  def add(self, engine):
    engine.fdtd.addmesh(x=self.bbox.pos.x, y=self.bbox.pos.y, z=self.bbox.pos.z, \
      x_span=self.bbox.size.x, y_span=self.bbox.size.y, z_span=self.bbox.size.z, \
      dx=(self.dx*1e6), dy=(self.dy*1e6), dz=(self.dz*1e6))
