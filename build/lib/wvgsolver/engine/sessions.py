from .base import Session
import os
import time
import threading
import glob
import logging
import tempfile
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

    self.fdtd = self.engine.lumapi.FDTD(hide=self.engine.hide)
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

  def set_sim_region(self, pos=None, size=None, boundaries={}):
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
