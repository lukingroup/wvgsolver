from ..geometry.sources import DipoleSource, ModeSource
from ..utils.linalg import Vec3, BBox
from ..utils.constants import F_EPSILON, AXIS_X, AXIS_Y, AXIS_Z, C_LIGHT
from ..analysis.procedures import FrequencySpectrum, SideWavePower, WaveEnergy, WaveProfile, Transmission
from .base import SimulationObject
from ..parse.plotables import Bandstructure, EField, Quasipotential
import logging
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import json

class UnitCell(SimulationObject):
  """This class represents a Nanophotonic unit cell as part of a cavity,
  and provides several simulations described below"""
  def __init__(self, structures=[], size=Vec3(1e-9), engine=None, load_path=None, metadata=None):
    """
    Parameters
    ----------
    structures : Structure or list
      The list of structures (or one structure) that make up this unit cell's geometry. The direction of
      periodicity of the unit cell is taken to be the x axis, so make sure you make the structures consistent
      with that
    size : Vec3
      The length of each side of the bounding box of this cell. The cell is assumed
      to be centered at the origin
    engine : Engine
      Engine to run simulations with
    load_path : str or None
      If provided, loads a unit cell from the given file path and populates all internal 
      data, including the list of structures, from that file
    metadata : any
      Any associated metadata
    """
    self._structures = structures if isinstance(structures, list) else [structures]
    self._size = size

    super().__init__(engine=engine, load_path=load_path, metadata=metadata)

    self._default_sim_type = "bandstructure"

  def _load_data(self, data):
    self._size = data["size"]
    self._structures = data["structures"]

  def _save_data(self):
    return {
      "size": self._size,
      "structures": self._structures
    }

  def get_structures(self, sim_type=None):
    if sim_type is None:
      return self._structures

    structures = []
    for s in self._structures:
      copy1_s = s.copy()
      copy2_s = s.copy()
      copy1_s.pos.x -= self._size.x / 2
      copy2_s.pos.x += self._size.x / 2
      structures.append(copy1_s)
      structures.append(copy2_s)

    return structures
  
  def get_size(self):
    return self._size

  def _simulate_bandgap(self, sess, freqs=(0.2e15, 0.6e15, 100000), min_gap=1e13, k=0.5, **kwargs):
    """Simulates and calculates the bandgap of this unit cell. It does this by simulating the 
    bandstructure at the given normalized wave vector (usually 0.5), finding the two highest peaks
    in the frequency domain around the region of interest, and returns them sorted.

    Parameters
    ----------
    sess : Session
      Used internally to perform the actual simulation
    freqs : tuple
      A tuple containing the starting frequency, ending frequency, and number of data points to perform
      the fourier transform over. This range should completely contain where you expect the bottom and top 
      of the band gap to be, as any bandgap returned will lie within it.
    min_gap : float
      Minimum bandgap that you're looking for. This ensures we don't accidently choose two points on the same 
      peak as the bandgap
    k : float
      Wave vector to calculate the bandgap at. Should usually be 0.5.
    **kwargs : dict
      Additional kwargs to pass to the bandstructure simulation call. See _simulate_bandstructure below
      for a description of the options
    
    Returns
    -------
    bandgap : tuple
      A two-element tuple, with the first element being the bottom of the bandgap, and the second being the top
    """
    res = self._simulate_bandstructure(sess, ks=(k, k, 1), freqs=freqs, **kwargs).data[0,:]
    freqs = np.linspace(*freqs)

    resolution = (np.max(freqs)-np.min(freqs))/len(freqs)
    distance = 15e12/resolution
    peaks, _ = find_peaks(res,prominence=2e4,distance=distance)
    #np.savetxt("freqs.txt",freqs)
    #np.savetxt("res.txt",res)

    #plt.plot(freqs,res)
    #plt.plot(freqs[peaks],res[peaks],'.')
    #plt.show()

    #first = freqs[np.argmax(res)]
    first = freqs[peaks[0]]
    second = freqs[peaks[1]]

    #bottom_ceil = max(0, int((first - min_gap - freqs[0]) * len(freqs)/(freqs[-1] - freqs[0])))
    #top_floor = min(len(freqs)-1, int((first + min_gap - freqs[0]) * len(freqs)/(freqs[-1] - freqs[0])))
    #if np.max(res[:bottom_ceil+1]) > np.max(res[top_floor:]):
    #  second = freqs[np.argmax(res[:bottom_ceil+1])]
    #else:
    #  second = freqs[np.argmax(res[top_floor:])+top_floor]

    #if first > second:
    #  return second, first
    return first, second

  def _simulate_bandstructure(self, sess, ks=(0, 0.5, 20), freqs=(0.2e15, 0.6e15, 100000), run_time=600e-15, \
    window_pos=0.5, TEonly=True, ndipoles=5, dipole_region=Vec3(1, 0, 0), dipole_directions=Vec3(0, 1, 0), \
    sim_size=3, analyze_region=0.1): 
    """Simulate the bandstructure of this unit cell by calculating the frequency domain left after an
    electric dipole pulse is allowed to dissipate through an infinite (implemented using boundary conditions)
    array of this unit cell, over a range of values for the wave vector.
    Specifically, the simulation consists of taking two of this unit cell side by side, and simulating a 
    region the size of the unit cell with Bloch boundary conditions on either end along the axis that the 
    two copies are stacked. A set of dipoles are randomly spread throughout the simulation region, and the
    simulation runs for the given time, collecting the frequency domain at the end.

    Parameters
    ----------
    sess : Session
      Used internally to perform the actual simulation
    ks : tuple
      The range of values of the wave vector to perform the simulation over, with the number of data points
      to break the range into as the third element. One simulation is run for each of the data points
    freqs : tuple
      The range of frequencies to perform the frequency analysis over, with the number of data points with
      which to perform the fourier transform as the third element
    run_time : float
      The run time of the simulation in seconds. The analysis is performed upon completion of this time
    window_pos : float
      A float between 0 and 1 indicating where to place the simulation window within the two stacked copies
      of the unit cell. A value of 0 means make the simulation window overlap exactly with the first copy, then
      as we increase it to 0.5 we move the window to be centered between the two copies, and then moving it to 1
      makes it overlap exactly with the second cell.
    TEonly : bool
      Whether or not to enforce TE only modes by applying an antisymmetric boundary condition along the y-axis.
    ndipoles : int
      The number of dipoles to use to generate the initial pulse
    dipole_region : Vec3 or float
      The size of the region within which the source dipoles are randomly placed, in relative units to the
      unit cell's size. For example, a value of Vec3(0.5, 0.2, 0.3) means the dipoles will be placed in a box sized
      0.5*x, 0.2*y, 0.3*z centered at the origin of the window, where x, y, and z are the dimensions of the unit cell.
      A single float value corresponds to a vector whose all three elements are equal to that value.
    dipole_directions : Vec3
      This indicates the relative weights of the three dimensions that the dipole axes can point along. Each dipole starts
      with a randomly chosen axis, and we multiply by dipole_directions to get the final axis. So if you wanted to force all
      dipoles to lie along the y-axis, you could provide Vec3(0, 1, 0)
    sim_size : Vec3 or float
      The size of the simulation region in the y and z dimensions (along the x dimension, it's equal to the x dimension of the unit
      cell because that's the dimension of periodicity), in units of the size of the unit cell. The x component of this vector is ignored, 
      so a value of Vec3(0, 1.5, 2) would correspond to a simulation region 1.5 times as large as the y-dimension of the unit cell
      in the y dimension, and 2 times as large as the z-dimension of the unit cell in the z dimension. A single float value 
      corresponds to a vector whose all three elements are equal to that value.
    analyze_region : Vec3 or float
      The size of the region within which frequency domain data will be gathered at a series of random points, in units
      of the size of the unit cell. A single float value corresponds to a vector whose all three elements are equal to that value.
    """
    ks = np.linspace(*ks)
    freqs = np.linspace(*freqs)
    
    sim_size = self._size*sim_size
    sim_size.x = self._size.x

    sess.set_sim_region(pos=Vec3(self._size.x * (window_pos - 0.5), 0, 0), size=sim_size, boundaries={
      "ymin": "antisymmetric" if TEonly else "pml"
    })
    sess.set_sim_time(run_time)

    nsweeps = len(ks)
    output = np.zeros((nsweeps, len(freqs)))
    for s in range(nsweeps):
      k = ks[s]
      sess.set_sim_region(boundaries={ "x": { "type": "bloch", "k": k } })
      dipoles = []
      for i in range(ndipoles):
        r = np.random.random_sample((3, )) - 0.5
        pos = Vec3(r[0] + (window_pos-0.5), 0.5*(r[1] + 0.5) if TEonly else r[1], r[2]) * self._size * dipole_region
        d = 2*(np.random.random_sample((3,)) - 0.5)
        dipoles.append(DipoleSource(frange=(freqs[0], freqs[-1]), pos=pos, axis=Vec3(d[0], d[1], d[2]) * dipole_directions, phase=(pos.x*k*2*np.pi/self._size.x)))
      sess.set_sources(dipoles)

      sweep = sess.run(FrequencySpectrum(BBox(
        Vec3(self._size.x * (window_pos - 0.5), self._size.y*0.5*(analyze_region.y if isinstance(analyze_region, Vec3) else analyze_region) if TEonly else 0, 0),
        self._size*analyze_region
      ), freqs))
      output[s,:] = sweep[:,0]
      logging.info("Sweep %d/%d completed (k=%f)" % (s + 1, nsweeps, k))
    
    return Bandstructure(output, (ks, freqs, self._size.x))

class Waveguide(SimulationObject):
  """This class represents a general waveguide structure"""
  def __init__(self, structures=[], size=Vec3(1e-6), engine=None, load_path=None, metadata=None):
    """
    Parameters
    ----------
    structures : list or Structure
      Structures that comprise this waveguide
    size : Vec3
      The size of the overal waveguide
    engine : Engine
      Engine to run simulations with
    load_path : str or None
      If provided, loads a waveguide from the given file path
    metadata : any
      Any associated metadata
    """
   
    self._structures = structures if isinstance(structures, list) else [structures] 
    self._size = size
    super().__init__(engine=engine, load_path=load_path, metadata=metadata)

    self._default_sim_type = "guidedness"

  def _save_data(self):
    return {
      "structures": self._structures
    }

  def _load_data(self, data):
    self._structures = data["structures"]

  def get_structures(self, sim_type=None):
    return self._structures
 
  def _getsize(self):
    return self._size

  def _simulate_guidedness(self, sess, target_freq=400e12, bboxes={}, simbbox=None, source_offset=0.5e-6, freq_span=0, freq_points=1, sim_time=200e-15, source_size=3, TEonly=True):
    """Simulates the guidedness of the cavity by using a mode source in the positive x direction and calculating 
    the transmission and reflection coefficients along the x axis. The sum of these coefficients indicates the 
    fraction of light from the mode source that remains guided along the x axis

    Parameters
    ----------
    sess : Session
      Used internally to perform the actual simulation
    bboxes : BBox or dict
      Specifies the bounding box(es) whose edges serve as the analysis regions for transmission calculations. Specifying
      a single BBox here will result in all 6 transmission monitors (2 for each axis) along the edges of that box. 
      Alternatively, you can specificy the boxes for specific axes by giving a dictionary with keys "x", "y", "z"
      and values being BBoxes. The default in any case is the simulation region.
    simbbox : BBox or None
      The explicit simulation bounding box. If None, use the waveguide size centered at 0.
    target_freq : float
      target frequency of the cavity simulation
    source_offset : float
      The x position offset of the mode source from the x min wall of bboxes["x"]
    freq_span : float
      The frequency span of the mode source
    sim_time : float
      The simulation time
    source_size : float or Vec3
      The size of the y and z dimensions of the mode source, in units of the size of the cavity in those dimensions
    freq_points : int
      The number of frequency points in the fourier transform of the transmission/reflection.
    TEonly : bool
      Whether or not to enforce TE only modes by applying an antisymmetric boundary condition along the y-axis.

    Returns
    -------
    transmission: dict
      Of the form { "x": float, "y": float, "z": float }.
      The sum of the transmission and reflection coefficients across the whole simulation time along the 3 axes
    """
    size = self._getsize()
    if simbbox is None:
      simbbox = BBox(Vec3(0), Vec3(1, 4, 4)*size)

    tbboxes = {
      "x": simbbox,
      "y": simbbox,
      "z": simbbox
    }
    if isinstance(bboxes, BBox):
      for a in tbboxes.keys():
        tbboxes[a] = bboxes
    elif isinstance(bboxes, dict):
      for a, b in enumerate(bboxes):
        tbboxes[a] = b

    sess.set_sim_region(size=simbbox.size, pos=simbbox.pos, boundaries={
      "ymin": "antisymmetric" if TEonly else "pml"
    })
    sess.set_sources(ModeSource(frange=(target_freq - freq_span * 0.5, target_freq + freq_span * 0.5),
      pos=Vec3(tbboxes["x"].pos.x - tbboxes["x"].size.x/2 + source_offset, 0, 0), size=(size * source_size)))
    sess.set_sim_time(sim_time)

    analysis = {
      "pxmin": Transmission(tbboxes["x"], AXIS_X, -1, target_freq, freq_span, freq_points),
      "pxmax": Transmission(tbboxes["x"], AXIS_X, 1, target_freq, freq_span, freq_points),
      "pymax": Transmission(tbboxes["y"], AXIS_Y, 1, target_freq, freq_span, freq_points),
      "pzmin": Transmission(tbboxes["z"], AXIS_Z, -1, target_freq, freq_span, freq_points),
      "pzmax": Transmission(tbboxes["z"], AXIS_Z, 1, target_freq, freq_span, freq_points),
    }
    if not TEonly:
      analysis["pymin"] = Transmission(tbboxes["y"], AXIS_Y, -1, target_freq, freq_span, freq_points),

    res = sess.run(analysis)

    return {
      "x": res["pxmax"] - res["pxmin"],
      "y": 2*res["pymax"] if TEonly else res["pymax"] - res["pymin"],
      "z": res["pzmax"] - res["pzmin"]
    }

  
class Cavity1D(Waveguide):
  """This class represents a 1D cavity, which consists of a series of unit cells stacked together, plus 
  additional structures (such as a beam). All simulations on this cavity are performed by stacking the list of
  unit cells in the x dimension, centering this stack around the origin, and adding any additional structures.
  """
  def __init__(self, unit_cells=[], structures=[], size=None, center_cell=None, center_shift=None, engine=None, load_path=None, metadata=None):
    """
    Parameters
    ----------
    unit_cells : list
      Array of unit cells in this cavity
    structures : list or Structure
      Additional structures that comprise this cavity
    size : Vec3
      If provided, this is taken to be the size of the cavity. If not provided, the size is calculated from the 
      stack of unit cells
    center_cell : int or None
      If provided, the center of the cavity will be taken to be the righthand (in the positive x direction) wall of 
      the unit cell at this index. If not provided, the center of the cavity will be the righthand wall of the 
      cell at index floor(n / 2) - 1 for n unit cells.
    center_shift : float or None
      An additional arbitrary value by which to shift the cavity when determining its center. A positive value 
      would shift the cavity to the right, making the center more towards the left of the unit cell stack. If not provided,
      this is taken to be 0 for an even number of cells, and -(center cell x size)/2 for an odd number of cells, which 
      ensures that with a value of center_cell=None the cavity center will be at the center of the stack.
    engine : Engine
      Engine to run simulations with
    load_path : str or None
      If provided, loads a cavity from the given file path and populates all internal 
      data, including the list of unit cells, from that file
    metadata : any
      Any associated metadata
    """
    
    self._unit_cells = unit_cells if isinstance(unit_cells, list) else [unit_cells]
    self._center_cell = center_cell
    self._center_shift = center_shift
    self._size_override = size

    super().__init__(structures=structures, engine=engine, load_path=load_path, size=size, metadata=metadata)

    self._default_sim_type = "resonance"
    self._no_sess_sims = ["quasipotential"]

  def get_unit_cells(self):
    return self._unit_cells

  def _load_data(self, data):
    self._structures = data["additional_structures"]
    unit_cells = []
    for s in data["unit_cells"]:
      c = UnitCell(engine=self.engine)
      c.load_data(s)
      unit_cells.append(c)
    self._unit_cells = unit_cells

  def _save_data(self):
    return {
      "additional_structures": self._structures,
      "unit_cells": [ c.save_data() for c in self._unit_cells ]
    }

  def get_structures(self, sim_type=None):
    if len(self._unit_cells) < 1:
      return [ s.copy() for s in self._structures ]
    
    center_cell = self._center_cell
    if center_cell is None:
      center_cell = int(len(self._unit_cells) / 2) - 1

    center_shift = self._center_shift
    if center_shift is None:
      center_shift = 0
      if len(self._unit_cells) % 2 == 1:
        center_shift = -self._unit_cells[int(len(self._unit_cells) / 2)].get_size().x / 2

    structs = []
    offset = Vec3(0)
    shift = 0
    for i in range(len(self._unit_cells)):
      c = self._unit_cells[i]
      xsize = c.get_size().x
      for s in c.get_structures():
        copy_s = s.copy()
        copy_s.pos = Vec3(offset) + copy_s.pos
        copy_s.pos.x += xsize / 2
        structs.append(copy_s)

      if i <= center_cell:
        shift += xsize
      offset.x += xsize

    for s in structs:
      s.pos.x -= shift - center_shift

    for s in self._structures:
      copy_s = s.copy()
      structs.append(copy_s)

    return structs

  def _getsize(self):
    if self._size_override:
      return self._size_override

    x = 0
    y = 0
    z = 0

    for c in self._unit_cells:
      size = c.get_size()
      x += size.x
      y = max(y, size.y)
      z = max(z, size.z)

    return Vec3(x, y, z)

  def _simulate_quasipotential(self, target_freq=400e12, **kwargs):
    """Simulates the quasipotential of the cavity by calculating the bandgap of each unit cell,
    and taking the minimum difference between the target frequency and each of the two
    edges of the gap. This distance can be negative, if the target frequency lies outside the gap.
    In this case, the quasipotential of that cell is negative.

    Parameters
    ----------
    target_freq : float
      target frequency of the cavity, this serves as 0 in the potential structure
    **kwargs : dict
      Additional arguments to pass to the unit cell's _simulate_bandgap call
    
    Returns
    -------
    Quasipotential
      A plotable Parser containing the quasipotential data. See parse/plotables for more info
    """
    output = []
    memo = []
    for cidx, c in enumerate(self._unit_cells):
      found = False
      for m in memo:
        if c.eq_structs(m[0]):
          output.append(m[1])
          found = True

      if not found:
        logging.info("Simulating cell %d/%d" % (cidx + 1, len(self._unit_cells)))
        gap = c.simulate("bandgap", **kwargs)
        r = target_freq - gap[0]
        output.append(r)
        memo.append((c, r))
    
    return Quasipotential(output)

  def _simulate_resonance(self, sess, target_freq=400e12, source_pulselength=60e-15, analyze_fspan=3e14, \
      analyze_time=590e-15, eref_time=80e-15, TEonly=True, sim_size=Vec3(2, 4, 4), energy_downsample=2):
    """Simulate the cavity's resonance frequency and obtain directional Q factors, mode volume, and 
    electric field profile data. The simulation comprises of one dipole polarized along the y-axis
    in the center of the cavity, and all simulation results are calculated in the last 20 fs 
    or so of the simulation. The electric dipole emits a pulse of a given length with a center
    frequency equal to the target frequency.

    Parameters
    ----------
    sess : Session
      Used internally to perform the actual simulation
    target_freq : float
      target frequency of the cavity simulation
    source_pulselength : float
      The length in seconds of the dipole pulse centered at the target frequency
    analyze_fspan : float
      The range of frequencies, centered at the target frequency, within which to look 
      for the resonance frequency
    analyze_time : float
      The time, in seconds, to simulation for before performing analysis. The total simulation time
      is a few dozens of femtoseconds greater than this
    eref_time : float
      The time, in seconds, at which to calculate the total energy in the cavity when determining
      the initial energy of the dipole pulse. This is then used to see how much energy is retained
      by the cavity.
    TEonly : bool
      Whether or not to enforce TE only modes by applying an antisymmetric boundary condition along the y-axis.
    sim_size : Vec3 or float
      The size of the simulation region in units of the size of the cavity. A single float value 
      corresponds to a vector whose all three elements are equal to that value.
    energy_downsample : int
      How much to downsample the electric field spatial profile data collection. Upping this value is handy
      to save memory when you have a very fine simulation mesh
    
    Returns
    -------
    dict
      A dict containing all the resulting data, of the form
      {
        "freq": float,
        "qxmin": float,
        "qxmax": float,
        "qymax": float,
        "qzmin": float,
        "qzmax": float,
        "vmode": float,
        "eremain": float,
        "xyprofile": Plotable EField object,
        "yzprofile": Plotable EField object
      }
      "freq" is the calculated resonance frequency, each "q..." result is a directional q factor,
      "vmode" is the volume of the cavity mode normalized to the cube of the resonance frequency,
      "xyprofile" and "yzprofile" are the electric field profiles in two planes that slice through 
      the center of the cavity, the first being along the x and y axes and the second along y and z.
      Finally, "eremain" is the approximate fraction of energy remaining in the simuation at the end.
    """
    size = self._getsize()
    analyze_frange = (target_freq - analyze_fspan / 2, target_freq + analyze_fspan / 2)

    bbox = BBox(Vec3(0), size * sim_size)

    sess.set_sim_region(size=bbox.size * Vec3(1.2, 1, 1), boundaries={
      "ymin": "antisymmetric" if TEonly else "pml"
    })
    sess.set_sources(DipoleSource(f=target_freq, pulse_length=source_pulselength, pulse_offset=2.1*source_pulselength, pos=Vec3(0), axis=Vec3(0, 1, 0)))
    sess.set_sim_time(analyze_time + 5 / target_freq)

    st = analyze_time

    spectrum_freqs = np.linspace(analyze_frange[0], analyze_frange[1], 100000)
    analysis = {
      "res": FrequencySpectrum(BBox(bbox.pos, size*Vec3(0.3, 0, 0)), spectrum_freqs),
      "eref": WaveEnergy(bbox, target_freq, eref_time, energy_downsample),
      "e": WaveEnergy(bbox, target_freq, st, energy_downsample),
      "pxmin": SideWavePower(bbox, AXIS_X, -1, st, target_freq),
      "pxmax": SideWavePower(bbox, AXIS_X, 1, st, target_freq),
      "pymax": SideWavePower(bbox, AXIS_Y, 1, st, target_freq),
      "pzmin": SideWavePower(bbox, AXIS_Z, -1, st, target_freq),
      "pzmax": SideWavePower(bbox, AXIS_Z, 1, st, target_freq),
      "xyprofile": WaveProfile(BBox(bbox.pos + Vec3(0, 0, 0.45*size.z), bbox.size), AXIS_Z, st, target_freq),
      "yzprofile": WaveProfile(bbox, AXIS_X, st, target_freq),
    }
    if not TEonly:
      analysis["pymin"] = SideWavePower(bbox, AXIS_Y, -1, st, target_freq)
    
    sess.run(analysis, analyze=False)

    res = sess.analyze("res")
    freq = spectrum_freqs[np.argmax(res)]

    e, edensity, index_x,_,_ = sess.analyze("e", freq = freq)
    eref,_,_,_,_ = sess.analyze("eref", freq = freq)

    vmode = (e / edensity) / (C_LIGHT / (freq * index_x))**3

    qxmin = 2*np.pi*freq*e/sess.analyze("pxmin")
    qxmax = 2*np.pi*freq*e/sess.analyze("pxmax")
    if not TEonly:
      qymin = 2*np.pi*freq*e/sess.analyze("pymin")
    qymax = 2*np.pi*freq*e/sess.analyze("pymax")
    qzmin = 2*np.pi*freq*e/sess.analyze("pzmin")
    qzmax = 2*np.pi*freq*e/sess.analyze("pzmax")
    xyprofile = sess.analyze("xyprofile")
    yzprofile = sess.analyze("yzprofile")

    output = {
      "freq": freq,
      "qxmin": qxmin,
      "qxmax": qxmax,
      "qymax": qymax,
      "qzmin": qzmin,
      "qzmax": qzmax,
      "vmode": vmode,
      "eremain": e,
      "xyprofile": EField(xyprofile, ("x", "y")),
      "yzprofile": EField(yzprofile, ("y", "z"))
    }
    if not TEonly:
      output["qymin"] = qymin
    return output
