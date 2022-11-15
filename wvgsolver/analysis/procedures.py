from .base import Analysis
from ..utils.misc import randstring
from ..utils.linalg import Vec3
from ..utils.constants import C_LIGHT, F_EPSILON, AXIS_X, AXIS_Y, AXIS_Z
import time
import numpy as np
try:
  import meep as mp
except:
  pass
from ..utils.meep import U_A, U_K, U_F, U_T

class FrequencySpectrum(Analysis):
  def __init__(self, bbox, freqs, nmonitors=5, start_time=0):
    self.nmonitors = nmonitors
    self.monitor_names = []
    for i in range(self.nmonitors):
      self.monitor_names.append(randstring())
    self.bbox = bbox
    self.freqs = freqs
    self.start_time = start_time

  def _setup_lumerical(self, sess):
    for i in range(self.nmonitors):
      r = np.random.random_sample((3, )) - 0.5
      pos = self.bbox.pos + Vec3(r[0], r[1], r[2]) * self.bbox.size
      sess.fdtd.addtime(name=self.monitor_names[i], monitor_type=1, x=pos.x, y=pos.y,
      z=pos.z, output_Hx=False, output_Hy=False, output_Hz=False, min_sampling_per_cycle=10000)

  def _cleanup_lumerical(self, sess):
    for i in range(self.nmonitors):
      sess.fdtd.select(self.monitor_names[i])
      sess.fdtd.delete()

  def _analyze_lumerical(self, sess):
    sess.fdtd.putv("f", self.freqs)
    script = """
    mnames = {%s};
    comps = {"Ex", "Ey", "Ez"};
    fd = 0;
    for(m=mnames) {
      t = pinch(getdata(m, 't'));
      for(comp=comps) {
        msignal = pinch(getdata(m, comp)) * exp(-0.5*(t-max(t)*0.5)^2/((max(t)*0.25)^2));
        fd = fd + abs(czt(msignal,t,2*pi*f)); 
      }
    }

    """ % (", ".join([ "'" + s + "'" for s in self.monitor_names ]))
    sess.fdtd.eval(script)

    return sess.fdtd.getv("fd")

  def _setup_eff1d(self, sess):
    pass
  
  def _cleanup_eff1d(self, sess):
    pass
  
  def _analyze_eff1d(self, sess):
    res = np.zeros(len(self.freqs))
    res[np.abs(self.freqs - sess.engine.mode_f).argmin()] = 1
    return res

class WaveEnergy(Analysis):
  def __init__(self, bbox, target_freq, start_time, downsample):
    self.bbox = bbox
    self.target_freq = target_freq
    self.start_time = start_time
    self.stop_time = start_time + 3 / target_freq
    self.downsample = downsample

    self.time_monitor = randstring()
    self.index_monitor = randstring()

  def _setup_lumerical(self, sess):
    time = sess.fdtd.addtime(name=self.time_monitor, monitor_type=8, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, output_Hx=False, output_Hy=False, output_Hz=False, down_sample_X=self.downsample,
      down_sample_Y=self.downsample, down_sample_Z=self.downsample)
    time.x_span = self.bbox.size.x
    time.y_span = self.bbox.size.y
    time.z_span = self.bbox.size.z
    time.stop_method = 2
    time.start_time = self.start_time
    time.stop_time = self.stop_time

    index = sess.fdtd.addindex(name=self.index_monitor, monitor_type=4, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, down_sample_X=self.downsample, down_sample_Y=self.downsample, down_sample_Z=self.downsample, use_source_limits=True)
    index.x_span = self.bbox.size.x
    index.y_span = self.bbox.size.y
    index.z_span = self.bbox.size.z

  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.time_monitor)
    sess.fdtd.delete()
    sess.fdtd.select(self.index_monitor)
    sess.fdtd.delete()

  def _analyze_lumerical(self, sess, freq):
    tm = self.time_monitor
    im = self.index_monitor

    sess.fdtd.putv("freq", freq)

    script = """
    im = "%s";
    index_x = pinch(getdata(im, "index_x"));
    index_y = pinch(getdata(im, "index_y"));
    index_z = pinch(getdata(im, "index_z"));
    
    tm = "%s";
    Ex = getdata(tm, "Ex");
    Ey = getdata(tm, "Ey");
    Ez = getdata(tm, "Ez");
    xv = getdata(tm, "x");
    yv = getdata(tm, "y");
    zv = getdata(tm, "z");
    t = getdata(tm, "t");

    Tpoints = length(t);

    E = matrix(Tpoints);
    E_maxdensity = matrix(Tpoints);
    for(i=1:Tpoints) {
        Edensity = (index_x^2 * (pinch(Ex, 4, i))^2 + index_y^2 * (pinch(Ey, 4, i))^2 + index_z^2 * (pinch(Ez, 4, i))^2);
        E(i) = integrate(Edensity, 1:3, xv, yv, zv);
        E_maxdensity(i) = max(Edensity);
    }

    E_avg = real(eps0 * 0.5 * (max(E) + min(E)));
    E_maxdensity_avg = real(eps0 * 0.5 * (max(E_maxdensity) + min(E_maxdensity)));
    index_x_max = max(index_x);
    index_y_max = max(index_y);
    index_z_max = max(index_z);
    """ % (im, tm)

    sess.fdtd.eval(script)
    return sess.fdtd.getv("E_avg"), sess.fdtd.getv("E_maxdensity_avg"), \
      sess.fdtd.getv("index_x_max"), sess.fdtd.getv("index_y_max"), \
      sess.fdtd.getv("index_z_max")
  
  def set_resobj(self, En):
    self.En = En
    #self.En = val.field_energy_in_box(
    #  center=mp.Vector3(z=self.bbox.pos.x/U_A),
    #  size=mp.Vector3(z=self.bbox.size.x/U_A)
    #)

  def _setup_eff1d(self, sess):
    self.region = {
      "time": self.start_time/U_T,
      "rtime": self.stop_time/U_T,
      "type": "energy",
      "args": (self.target_freq/U_F, 0, 1, mp.EnergyRegion(
        center=mp.Vector3(z=self.bbox.pos.x/U_A),
        size=mp.Vector3(z=self.bbox.size.x/U_A)
      )),
      "kwargs": {},
      "callback": self.set_resobj
    }
    sess.add_analyze_region(self.region)
    #sess.add_analyze_func({
    #  "time": self.start_time/U_T,
    #  "func": self.set_resobj
    #})
  
  def _cleanup_eff1d(self, sess):
    sess.remove_analyze_region(self.region)
    del self.region
  
  def _analyze_eff1d(self, sess, freq):
    e = mp.get_total_energy(self.En)[0]
    del self.En
    return e, np.inf, 1, 1, 1

class WaveProfile(Analysis):
  def __init__(self, bbox, axis, start_time, target_freq):
    self.bbox = bbox
    self.axis = axis
    self.start_time = start_time
    self.target_freq = target_freq
    self.time_monitor = randstring()
    self.index_monitor = randstring()
  
  def _setup_lumerical(self, sess):
    index_map = {}
    index_map[AXIS_X] = 1
    index_map[AXIS_Y] = 2
    index_map[AXIS_Z] = 3

    time_map = {}
    time_map[AXIS_X] = 5
    time_map[AXIS_Y] = 6
    time_map[AXIS_Z] = 7

    time = sess.fdtd.addtime(name=self.time_monitor, monitor_type=time_map[self.axis], x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, output_Hx=False, output_Hy=False, output_Hz=False)
    time.stop_method = 2
    time.start_time = self.start_time
    # Allow for one period of the resonance oscillation
    time.stop_time = self.start_time + 1 / self.target_freq

    index = sess.fdtd.addindex(name=self.index_monitor, monitor_type=index_map[self.axis], x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, use_source_limits=True)

    if self.axis == AXIS_X:
      time.y_span = self.bbox.size.y
      time.z_span = self.bbox.size.z
      index.y_span = self.bbox.size.y
      index.z_span = self.bbox.size.z
    elif self.axis == AXIS_Y:
      time.x_span = self.bbox.size.x
      time.z_span = self.bbox.size.z
      index.x_span = self.bbox.size.x
      index.z_span = self.bbox.size.z
    else:
      time.x_span = self.bbox.size.x
      time.y_span = self.bbox.size.y
      index.x_span = self.bbox.size.x
      index.y_span = self.bbox.size.y
  
  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.time_monitor)
    sess.fdtd.delete()
    sess.fdtd.select(self.index_monitor)
    sess.fdtd.delete()
  
  def _analyze_lumerical(self, sess):
    tm = self.time_monitor
    im = self.index_monitor

    d_map = {}
    d_map[AXIS_X] = ["y", "z"]
    d_map[AXIS_Y] = ["x", "z"]
    d_map[AXIS_Z] = ["x", "y"]

    script = """
    im = "%s";
    index_x = pinch(getdata(im, "index_x"));
    index_y = pinch(getdata(im, "index_y"));
    index_z = pinch(getdata(im, "index_z"));
    
    tm = "%s";
    Ex = pinch(getdata(tm, "Ex"));
    Ey = pinch(getdata(tm, "Ey"));
    Ez = pinch(getdata(tm, "Ez"));
    d1 = pinch(getdata(tm, "%s"));
    d2 = pinch(getdata(tm, "%s"));
    t = pinch(getdata(tm, "t"));

    Tpoints = length(t);

    EdensityX = matrix(length(d1), length(d2));
    EdensityY = matrix(length(d1), length(d2));
    EdensityZ = matrix(length(d1), length(d2));
    for(i=1:Tpoints) {
        EdensityX = EdensityX + index_x^2 * (pinch(Ex, 3, i))^2;
        EdensityY = EdensityY + index_y^2 * (pinch(Ey, 3, i))^2;
        EdensityZ = EdensityZ + index_z^2 * (pinch(Ez, 3, i))^2;
    }
    EdensityX = real(eps0 * EdensityX / Tpoints);
    EdensityY = real(eps0 * EdensityY / Tpoints);
    EdensityZ = real(eps0 * EdensityZ / Tpoints);
    """ % (im, tm, d_map[self.axis][0], d_map[self.axis][1])

    sess.fdtd.eval(script)
    return sess.fdtd.getv("EdensityX"), sess.fdtd.getv("EdensityY"), \
      sess.fdtd.getv("EdensityZ"), np.squeeze(sess.fdtd.getv("d1")), \
      np.squeeze(sess.fdtd.getv("d2")), np.real(sess.fdtd.getv("index_x"))
  
  def _setup_eff1d(self, sess):
    pass
  
  def _cleanup_eff1d(self, sess):
    pass
  
  def _analyze_eff1d(self, sess):
    return None, None, None, None, None, None

class SideWavePower(Analysis):
  def __init__(self, bbox, axis, side, start_time, target_freq):
    self.bbox = bbox
    self.axis = axis
    self.side = side
    self.start_time = start_time
    self.stop_time = start_time + 3/target_freq
    self.target_freq = target_freq
    self.monitor_name = randstring()

  def _setup_lumerical(self, sess):
    type_map = {}
    type_map[AXIS_X] = 5
    type_map[AXIS_Y] = 6
    type_map[AXIS_Z] = 7

    time = sess.fdtd.addtime(name=self.monitor_name, monitor_type=type_map[self.axis], x=self.bbox.pos.x,
       y=self.bbox.pos.y, z=self.bbox.pos.z, output_Ex=False, output_Ey=False, output_Hx=False,
       output_Hy=False, output_Hz=False, output_Ez=False, output_power=True)
    time.stop_method = 2
    time.start_time = self.start_time
    time.stop_time = self.stop_time 
    
    if self.axis == AXIS_X:
      time.x = time.x + self.side * self.bbox.size.x / 2
      time.y_span = self.bbox.size.y
      time.z_span = self.bbox.size.z
    elif self.axis == AXIS_Y:
      time.x_span = self.bbox.size.x
      time.y = time.y + self.side * self.bbox.size.y / 2
      time.z_span = self.bbox.size.z
    else:
      time.x_span = self.bbox.size.x
      time.y_span = self.bbox.size.y
      time.z = time.z + self.side * self.bbox.size.z / 2

  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.monitor_name)
    sess.fdtd.delete()

  def _analyze_lumerical(self, sess):
    data = np.real(sess.fdtd.getdata(self.monitor_name, "power"))
    return self.side * 0.5 * (np.max(data) + np.min(data))
  
  def set_resobj(self, F):
    self.F = F

  def _setup_eff1d(self, sess):
    if self.axis != AXIS_X:
      return

    self.region = {
      "time": self.start_time/U_T,
      "rtime": self.stop_time/U_T, 
      "type": "flux",
      "args": (self.target_freq/U_F, 0, 1, mp.FluxRegion(
        center=mp.Vector3(z=(self.bbox.pos.x + self.side*self.bbox.size.x/2)/U_A),
        size=mp.Vector3()
      )),
      "kwargs": {},
      "callback": self.set_resobj
    }
    sess.add_analyze_region(self.region)
  
  def _cleanup_eff1d(self, sess):
    if hasattr(self, "region"):
      sess.remove_analyze_region(self.region)
      del self.region
  
  def _analyze_eff1d(self, sess):
    if not hasattr(self, "F"):
      return np.inf
    f = mp.get_fluxes(self.F)[0]
    del self.F
    return U_F*f*self.side

class Transmission(Analysis):
  def __init__(self, bbox, axis, side, target_freq, freq_span, freq_points=20, return_spectrum=False):
    self.bbox = bbox
    self.axis = axis
    self.side = side
    self.target_freq = target_freq
    self.freq_span = freq_span
    self.freq_points = freq_points
    self.return_spectrum = return_spectrum
    self.monitor_name = randstring()

  def _setup_lumerical(self, sess):
    type_map = {}
    type_map[AXIS_X] = 5
    type_map[AXIS_Y] = 6
    type_map[AXIS_Z] = 7

    power = sess.fdtd.addpower(name=self.monitor_name, monitor_type=type_map[self.axis], x=self.bbox.pos.x,
       y=self.bbox.pos.y, z=self.bbox.pos.z, output_Ex=False, output_Ey=False, output_Hx=False,
       output_Hy=False, output_Hz=False, output_Ez=False, output_power=True, output_Px=False, output_Py=False,
       output_Pz=False, override_global_monitor_settings=True)
    power.use_source_limits = False
    power.wavelength_center = C_LIGHT / self.target_freq
    power.wavelength_span = C_LIGHT / (self.target_freq - 0.5*self.freq_span) - C_LIGHT / (self.target_freq + 0.5*self.freq_span)
    power.frequency_points = self.freq_points
    
    if self.axis == AXIS_X:
      power.x = power.x + self.side * self.bbox.size.x / 2
      power.y_span = self.bbox.size.y
      power.z_span = self.bbox.size.z
    elif self.axis == AXIS_Y:
      power.x_span = self.bbox.size.x
      power.y = power.y + self.side * self.bbox.size.y / 2
      power.z_span = self.bbox.size.z
    else:
      power.x_span = self.bbox.size.x
      power.y_span = self.bbox.size.y
      power.z = power.z + self.side * self.bbox.size.z / 2

  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.monitor_name)
    sess.fdtd.delete()

  def _analyze_lumerical(self, sess):
    sess.fdtd.putv("freq", self.target_freq)
    
    script = """
    pm = "%s";
    f = pinch(getdata(pm, "f"));
    T = transmission(pm);
    Tf = T(find(f, freq));
    """ % self.monitor_name
    sess.fdtd.eval(script)

    if self.return_spectrum:
      return {
        "f": sess.fdtd.getv("f"),
        "T": sess.fdtd.getv("T")
      }
    else:
      return sess.fdtd.getv("Tf")
  
  def _setup_eff1d(self, sess):
    if self.axis != AXIS_X:
      return

    self.regions = [{
      "time": 0,
      "rtime": -1, 
      "type": "flux",
      "args": (self.target_freq/U_F, 0, 1, mp.FluxRegion(
        center=mp.Vector3(z=(self.bbox.pos.x - self.bbox.size.x/2)/U_A),
        size=mp.Vector3()
      )),
      "kwargs": {},
      "callback": self.set_resobj1
    },{
      "time": 0,
      "rtime": -1, 
      "type": "flux",
      "args": (self.target_freq/U_F, 0, 1, mp.FluxRegion(
        center=mp.Vector3(z=(self.bbox.pos.x + self.bbox.size.x/2)/U_A),
        size=mp.Vector3()
      )),
      "kwargs": {},
      "callback": self.set_resobj2
    }]
    for r in self.regions:
      sess.add_analyze_region(r)

  def set_resobj1(self, F):
    self.F1 = F
  
  def set_resobj2(self, F):
    self.F2 = F
  
  def _cleanup_eff1d(self, sess):
    if hasattr(self, "regions"):
      for r in self.regions:
        sess.remove_analyze_region(r)
      del self.regions
  
  def _analyze_eff1d(self, sess):
    if not hasattr(self, "F1") or not hasattr(self, "F2"):
      return 0
    f1 = mp.get_fluxes(self.F1)[0]
    f2 = mp.get_fluxes(self.F2)[0]
    print(f1, f2)
    del self.F1, self.F2
    if self.side == 1:
      return U_F*f2 / (f2 - f1/2)
    return -U_F*(f1/2) / (f2 - f1/2)

class Fields(Analysis):
  def __init__(self, bbox, ndims=2, axis=AXIS_Z, start_time=0, downsample=1):
    self.monitor_name = randstring()
    self.bbox = bbox
    self.start_time = start_time
    self.ndims = ndims
    self.axis = axis
    self.downsample = downsample

  def _setup_lumerical(self, sess):
    type_map = {}
    type_map[AXIS_X] = 1
    type_map[AXIS_Y] = 2
    type_map[AXIS_Z] = 3

    mtype = 1
    if self.ndims == 1:
      mtype = type_map[self.axis] + 1
    elif self.ndims == 2:
      mtype = type_map[self.axis] + 4
    elif self.ndims == 3:
      mtype = 8

    time = sess.fdtd.addtime(name=self.monitor_name, monitor_type=mtype, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, output_Hx=False, output_Hy=False, output_Hz=False, output_Px=False, output_Py=False, output_Pz=False)
    time.stop_method = 1
    time.start_time = self.start_time
   
    if self.ndims == 1:
      if self.axis == AXIS_X:
        time.x_span = self.bbox.size.x
        time.down_sample_X = self.downsample
      elif self.axis == AXIS_Y:
        time.y_span = self.bbox.size.y
        time.down_sample_Y = self.downsample
      else:
        time.z_span = self.bbox.size.z
        time.down_sample_Z = self.downsample
    elif self.ndims == 2:
      if self.axis == AXIS_X:
        time.x = self.bbox.pos.x + self.bbox.size.x / 2
        time.y_span = self.bbox.size.y
        time.z_span = self.bbox.size.z
        time.down_sample_Y = self.downsample
        time.down_sample_Z = self.downsample
      elif self.axis == AXIS_Y:
        time.x_span = self.bbox.size.x
        time.y = self.bbox.pos.y + self.bbox.size.y / 2
        time.z_span = self.bbox.size.z
        time.down_sample_X = self.downsample
        time.down_sample_Z = self.downsample
      else:
        time.x_span = self.bbox.size.x
        time.y_span = self.bbox.size.y
        time.z = self.bbox.pos.z + self.bbox.size.z / 2
        time.down_sample_X = self.downsample
        time.down_sample_Y = self.downsample
    elif self.ndims == 3:
      time.x_span = self.bbox.size.x
      time.y_span = self.bbox.size.y
      time.z_span = self.bbox.size.z
      time.down_sample_X = self.downsample
      time.down_sample_Y = self.downsample
      time.down_sample_Z = self.downsample

  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.monitor_name)
    sess.fdtd.delete()

  def _analyze_lumerical(self, sess):
    fields = ["Ex", "Ey", "Ez"]

    return np.stack(
      list(filter(
        lambda a: isinstance(a, np.ndarray) and a.size != 0,
        [ np.real(sess.fdtd.getdata(self.monitor_name, f)) for f in fields ]
      ))
    )
  
  def _setup_eff1d(self, sess):
    pass
  
  def _cleanup_eff1d(self, sess):
    pass
  
  def _analyze_eff1d(self, sess):
    pass

class Index(Analysis):
  def __init__(self, bbox, ndims=3, norm_axis=AXIS_Z, axis=AXIS_X, downsample=1):
    self.monitor_name = randstring()
    self.bbox = bbox
    self.ndims = ndims
    self.norm_axis= norm_axis
    self.axis = axis
    self.downsample = downsample

  def _setup_lumerical(self, sess):
    type_map = {}
    type_map[AXIS_X] = 1
    type_map[AXIS_Y] = 2
    type_map[AXIS_Z] = 3

    mtype = 4
    if self.ndims == 2:
      mtype = type_map[self.norm_axis]

    index = sess.fdtd.addindex(name=self.monitor_name, monitor_type=mtype, x=self.bbox.pos.x, y=self.bbox.pos.y, z=self.bbox.pos.z)
   
    if self.ndims == 2:
      if self.norm_axis == AXIS_X:
        index.x = self.bbox.pos.x + self.bbox.size.x / 2
        index.y_span = self.bbox.size.y
        index.z_span = self.bbox.size.z
        time.down_sample_Y = self.downsample
        time.down_sample_Z = self.downsample
      elif self.norm_axis == AXIS_Y:
        index.x_span = self.bbox.size.x
        index.y = self.bbox.pos.y + self.bbox.size.y / 2
        index.z_span = self.bbox.size.z
        time.down_sample_X = self.downsample
        time.down_sample_Z = self.downsample
      else:
        index.x_span = self.bbox.size.x
        index.y_span = self.bbox.size.y
        index.z = self.bbox.pos.z + self.bbox.size.z / 2
        time.down_sample_X = self.downsample
        time.down_sample_Y = self.downsample
    else:
      index.x_span = self.bbox.size.x
      index.y_span = self.bbox.size.y
      index.z_span = self.bbox.size.z
      time.down_sample_X = self.downsample
      time.down_sample_Y = self.downsample
      time.down_sample_Z = self.downsample

  def _cleanup_lumerical(self, sess):
    sess.fdtd.select(self.monitor_name)
    sess.fdtd.delete()

  def _analyze_lumerical(self, sess):
    vmap = {}
    vmap[AXIS_X] = "x"
    vmap[AXIS_Y] = "y"
    vmap[AXIS_Z] = "z"

    return sess.fdtd.getdata(self.monitor_name, "index_" + vmap[self.axis])
  
  def _setup_eff1d(self, sess):
    pass
  
  def _cleanup_eff1d(self, sess):
    pass
  
  def _analyze_eff1d(self, sess):
    pass
