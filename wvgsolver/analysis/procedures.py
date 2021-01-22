from .base import Analysis
from ..utils.misc import randstring
from ..utils.constants import EPS0, AXIS_X, AXIS_Y, AXIS_Z
import time
import numpy as np

class FindResonances(Analysis):
  def __init__(self, bbox, frange, start_time=0):
    self.monitor_name = randstring()
    self.bbox = bbox
    self.frange = frange
    self.start_time = start_time

  def _setup_lumerical(self, sess):
    time = sess.fdtd.addtime(name=self.monitor_name, monitor_type=1, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, output_Ex=False, output_Ey=False, output_Hx=False, output_Hy=False,
      output_Hz=False, start_time=self.start_time)

  def _analyze_lumerical(self, sess):
    mn = self.monitor_name
    script = """
    signal = pinch(getdata('%s', 'Ez'));
    t = pinch(getdata('%s', 't'));

    res = findresonances(t, signal, [%f, %f]);
    """ % (mn, mn, self.frange[0], self.frange[1])

    sess.fdtd.eval(script)
    return sess.fdtd.getv("res")

class WaveEnergy(Analysis):
  def __init__(self, bbox, target_freq, start_time):
    self.bbox = bbox
    self.target_freq = target_freq
    self.start_time = start_time

    self.time_monitor = randstring()
    self.index_monitor = randstring()

  def _setup_lumerical(self, sess):
    time = sess.fdtd.addtime(name=self.time_monitor, monitor_type=8, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z, output_Hx=False, output_Hy=False, output_Hz=False)
    time.x_span = self.bbox.size.x
    time.y_span = self.bbox.size.y
    time.z_span = self.bbox.size.z
    time.stop_method = 2
    time.start_time = self.start_time
    # Allow for one wavelength
    time.stop_time = self.start_time + 1 / self.target_freq

    index = sess.fdtd.addindex(name=self.index_monitor, monitor_type=4, x=self.bbox.pos.x, y=self.bbox.pos.y,
      z=self.bbox.pos.z)
    index.x_span = self.bbox.size.x
    index.y_span = self.bbox.size.y
    index.z_span = self.bbox.size.z
    index.use_source_limits = False
    index.wavelength_center = 1 / self.target_freq
    index.wavelength_span = 0

  def _analyze_lumerical(self, sess, freq):
    tm = self.time_monitor
    im = self.index_monitor

    script = """
    index_x = pinch(getdata("%s", "index_x"));
    index_y = pinch(getdata("%s", "index_y"));
    index_z = pinch(getdata("%s", "index_z"));
    
    Ex = getdata("%s", "Ex");
    Ey = getdata("%s", "Ey");
    Ez = getdata("%s", "Ez");
    xv = getdata("%s", "x");
    yv = getdata("%s", "y");
    zv = getdata("%s", "z");
    t = getdata("%s", "t");

    dt = t(2) - t(1);
    Tpoints = round(0.5 * %f / dt) + 1;
    tv = t(1:Tpoints);

    E_avg = 0.0;
    for(i=1:length(tv)) {
        Edensity = (index_x^2 * (pinch(Ex, 4, i))^2 + index_y^2 * (pinch(Ey, 4, i))^2 + index_z^2 * (pinch(Ez, 4, i))^2);
        E_avg = E_avg + integrate(Edensity, 1:3, xv, yv, zv);
    }

    E_avg = real(%s * E_avg / Tpoints);
    """ % (im, im, im, tm, tm, tm, tm, tm, tm, tm, 1 / freq, "{:e}".format(EPS0))

    sess.fdtd.eval(script)
    return sess.fdtd.getv("E_avg")

class SidePower(Analysis):
  def __init__(self, bbox, axis, side, start_time, end_time):
    self.bbox = bbox
    self.axis = axis
    self.side = side
    self.start_time = start_time
    self.end_time = end_time
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
    time.stop_time = self.end_time
    
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

  def _analyze_lumerical(self, sess):
    return self.side * np.real(sess.fdtd.getdata(self.monitor_name, "power"))
