import uuid
from abc import ABC, abstractmethod
from ..utils.misc import hasmethod
from ..engine import getDefaultEngine
from collections import Counter
import pickle
import copy
import datetime

class SimulationObject(ABC):
  """
  Anything for which you want to run a simulation on is a SimulationObject.
  All SimulationObjects have several things in common:
  - Are primarily defined by a list of Structures, which represent the 3D geometry of 
    this object that is then used in simulation
  - Have the ability to run simulations using an Engine 
  - Have the ability to save and load themselves to and from a file
  - Have automatic storing and cataloguing of simulation results

  SimulationObjects follow a simple paradigm - once you've created them, don't modify their
  parameters (structures, unit cells, etc). This ensures that all the simulation results associated with
  a particular instance of a SimulationObject come from simulating the same object.
  """
  def __init__(self, engine=None, load_path=None, metadata=None):
    """
    Parameters
    ----------
    engine : Engine or None
      The engine to run the simulations of this object on. If not provided, we default to an
      instance of LumericalEngine
    load_path : str or None
      The file path to load this object from. If this parameter is not None, then the object from
      that file is loaded and you are not expected to modify the object further
    metadata : any
      Any associated metadata
    """
    if engine is None:
      engine = getDefaultEngine()

    self.engine = engine
    self._simulate_results = {}
    self._default_sim_type = ""
    self._no_sess_sims = []
    self.save_path = load_path
    self.metadata = metadata

    if load_path is not None:
      self.load(load_path)

  @abstractmethod
  def _save_data(self):
    pass
  
  @abstractmethod
  def _load_data(self, data):
    pass

  def save_data(self):
    """Get all data associated with this object

    Returns
    -------
    data : dict
      A dictionary of the form
      {
        ... Object-specific parameters
        "simulate_results": {
          "simulation_type": [
            {
              "res": dict of results (this is the main set of simulation results),
              "sess_res": A session-specific result object (this points to the Lumerical FSP file
                used for this simulation with the Lumerical Engine),
              "kwargs": Keyword arguments used for this simulation call,
              "engine": A string description of the engine used for this call,
              "start_t": Start time of this call,
              "end_t": End time of this call,
              "status": "succeeded" or "failed"
            }
          ]
        }
      }
    """
    data = self._save_data()
    data["simulate_results"] = {}
    data["metadata"] = copy.copy(self.metadata)
    for t in self._simulate_results:
      data["simulate_results"][t] = []
      for r in self._simulate_results[t]:
        data["simulate_results"][t].append(copy.copy(r))
        if "sess_res" in r:
          data["simulate_results"][t][-1]["sess_res"] = r["sess_res"].get_clean()
    return data
  
  def load_data(self, data):
    """Takes a set of data, usually returned by a call to save_data on the same type of object,
    and replaces all internal data in this object with the new data

    Parameters
    ----------
    data : dict
      Same format as the return value of save_data
    """
    self._load_data(data)
    self._simulate_results = data["simulate_results"]
    self.metadata = data["metadata"] if "metadata" in data else None
    for t in self._simulate_results:
      for r in self._simulate_results[t]:
        if "sess_res" in r:
          r["sess_res"].set_engine(self.engine)

  def save(self, fpath=None):
    """Save this object to a file path. The object then keeps track of this file path and re-saves
    itself to that file after every simulation

    Parameters
    ----------
    fpath : str
      The file path to save to. If None, attempts to use the previously set file path
    """
    if fpath is not None:
      self.save_path = fpath
    elif self.save_path is None:
      raise ValueError("No save path set, please provide one")
    pickle.dump(self.save_data(), open(self.save_path, "wb"))
  
  def load(self, fpath):
    """Load an object from a file path and replace all internal data with that from the given file.
    This does not keep track of the file path for future saves

    Parameters
    ----------
    fpath : str
    """
    self.load_data(pickle.load(open(fpath, "rb")))

  @abstractmethod
  def get_structures(self, sim_type=None):
    """Get the full list of structures in this object"""
    return []
  
  def _check_sim_type(self, t): 
    sim_type = "_simulate_" + t
    if not hasmethod(self, sim_type):
      raise NotImplementedError("simulation type '%s' not found for object '%s'" % \
        (t, type(self).__name__))
    return sim_type

  def get_results(self, sim_type=None, started_after=None, started_before=None, 
      ended_after=None, ended_before=None, status=None, **kwargs):
    """Get a filtered list of the results of past simulations on this object. When called with no arguments,
    this returns all results from all simulations on this object

    Parameters
    ----------
    sim_type : str
      If provided, only returns results of simulations of this type
    started_after : datetime.datetime
      Only return results of simulations that ran after this time
    started_before : datetime.datetime
      Only return results of simulations that started before this time
    ended_after : datetime.datetime
      Only return results of simulations that ended after this time
    ended_before : datetime.datetime
      Only return results of simulations that ended before this time
    status : {"succeeded", "failed", None}
      Only return results that ended with this status
    **kwargs : dict
      Only return results of simulations that we ran with a superset of the kwargs provided

    Returns
    -------
    results : dict or list
      If sim_type is None, then a dict of the form
      {
        "simulation_type": [
            {
              "res": dict of results (this is the main set of simulation results),
              "sess_res": A session-specific result object (this points to the Lumerical FSP file
                used for this simulation with the Lumerical Engine),
              "kwargs": Keyword arguments used for this simulation call,
              "engine": A string description of the engine used for this call,
              "start_t": Start time of this call,
              "end_t": End time of this call,
              "status": "succeeded" or "failed"
            }
          ]
        }
      }
      is returned. If sim_type is given, then a list of the same form as the "simulation_type"
      element is returned, for the simulation type corresponding to sim_type
    """
    if sim_type is None:
      return self._simulate_results
  
    if not sim_type in self._simulate_results:
      return []
    if not kwargs:
      res = self._simulate_results[sim_type]
    else:
      res =[
        results for results in self._simulate_results[sim_type] if 
          len(set(kwargs.items()) - set(results["kwargs"].items())) == 0
      ]
    
    if started_after is not None:
      res = filter(lambda r: r["start_t"] >= started_after, res)
    if started_before is not None:
      res = filter(lambda r: r["start_t"] <= started_before, res)
    
    if ended_after is not None:
      res = filter(lambda r: r["end_t"] >= ended_after, res)
    if ended_before is not None:
      res = filter(lambda r: r["end_t"] <= ended_before, res)
    
    if status is not None:
      res = filter(lambda r: r["status"] == status, res)

    return list(res)

  def add_simulation(self, t, func):
    """Adds a custom simulation.

    Parameters
    ----------
    t : str
      The simulation type being added. Can by any string representing a short 1-word or so 
      description of the simuation
    func : callable
      The simulation function
    """

    setattr(self, "_simulate_" + t, lambda *args, **kwargs: func(self, *args, **kwargs))
      

  def simulate(self, t=None, save=True, mesh_regions=[], **kwargs):
    """This is the main function that is called to run a simulation.
    Note: The values of the additional keyword arguments here can only be 
    Python primitives or tuples/lists/dicts of primitives, or Vec3's

    Parameters
    ----------
    t : str or None
      The simulation type to run. Must be implemented by the subclass of SimulationObject that the
      current instance belongs to. For example, if you want to run a simulation of type "resonance", 
      then you must implement _simulate_resonance in the object subclass. If None, then we use
      the default simulation type, which is also specified by the subclass
    save : boolean
      Whether or not to save the object automatically after running the simulation. This uses
      the last fpath passed to a call on save(), or if there was no such call previously, it 
      doesn't save the object at all
    **kwargs : dict
      Keyword arguments passed to the simulation implementation. Check simulation/objects.py
      for a description of the options for each simulation type
    """
    if t is None:
      t = self._default_sim_type

    sim_type = self._check_sim_type(t)

    args = []
    err = None
    sess_res = None
    session = None
    start_time = datetime.datetime.now()
    try:
      if t not in self._no_sess_sims:
        session = self.engine.new_session()
        session.set_structures(self.get_structures(t))
        session.set_mesh_regions(mesh_regions)

        args.append(session)

      res = getattr(self, sim_type)(*args, **kwargs)
    except Exception as e:
      err = e
    end_time = datetime.datetime.now()

    if t not in self._no_sess_sims and session:
      sess_res = session.get_postrunres()
      session.close()

    if t not in self._simulate_results:
      self._simulate_results[t] = []
    self._simulate_results[t].append({
      "kwargs": kwargs,
      "engine": str(self.engine),
      "start_t": start_time,
      "end_t": end_time,
      "status": "succeeded" if err is None else "failed"
    })
    if err is None:
      self._simulate_results[t][-1]["res"] = res
    if sess_res:
      self._simulate_results[t][-1]["sess_res"] = sess_res

    if save and self.save_path:
      self.save(self.save_path)

    if err is not None:
      raise err

    return res

  def eq_structs(self, other, sim_type=None):
    return Counter(self.get_structures(sim_type)) == Counter(other.get_structures(sim_type))
