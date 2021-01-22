from .base import Session
import os
import time
import threading
import glob
import logging

def lumericalLogFile(dir_path, name, stop):
  fps = {}
  while not stop.is_set():
    for f in glob.iglob(os.path.join(dir_path, name + "*.log")):
      try:
        base_name = os.path.basename(f)
        log_name = base_name[len(name):base_name.index(".log")]
        if log_name[0] == "_":
          log_name = log_name[1:]

        first = False
        if f not in fps:
          fps[f] = open(f, "r")
#          first = True

        while True:
          line = fps[f].readline()
          if line:
            if not first:
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

    self.structures = "structures"
    self.sources = "sources"
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

    self.working_path = os.path.join(self.engine.working_path, self.name)

    if not os.path.isdir(self.working_path):
      os.mkdir(self.working_path)

    self.fdtd = self.engine.lumapi.FDTD()
    self.fdtd.addstructuregroup(name=self.structures)
    self.fdtd.addgroup(name=self.sources)
    self.sim_region = self.fdtd.addfdtd(mesh_accuracy=4)

  def close(self):
    if self.fdtd is not None:
      self.fdtd.close()

    self.fdtd = None
    self.sim_region = None
  
  def _set_structures(self, structs=[]):
    self.fdtd.switchtolayout()
    self.fdtd.groupscope(self.structures)
    self.fdtd.deleteall()

    for s in structs:
      s.add(self)
    
    self.fdtd.groupscope("::model")
  
  def _set_sources(self, sources=[]):
    self.fdtd.switchtolayout()
    self.fdtd.groupscope(self.sources)
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
      else:
        mapped_key = self.boundary_keys_map[key]
        if not self.fdtd.ispropertyactive(self.sim_region.name, mapped_key):
          self.sim_region[mapped_key[0] + " min bc"] = 1
        self.sim_region[mapped_key] = self.boundary_values_map[val]
  
  def _prerun(self):
    self.fdtd.switchtolayout()
          
  def _runsim(self, options={}):
    self.fdtd.switchtolayout()
    self.fdtd.save(os.path.join(self.working_path, self.name + ".fsp"))

    stop_logging = threading.Event()
    if not "silent" in options or not options["silent"]:
      threading.Thread(target=lumericalLogFile, args=(self.working_path, self.name, stop_logging)).start()

    self.fdtd.run()

    stop_logging.set()
