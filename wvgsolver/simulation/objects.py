from ..geometry.sources import DipoleSource
from ..utils.linalg import Vec3, BBox
from ..utils.constants import F_EPSILON, AXIS_X, AXIS_Y, AXIS_Z
from ..analysis.procedures import FindResonances, SidePower, WaveEnergy
from .base import SimulationObject
import logging
import numpy as np

class UnitCell(SimulationObject):
  def __init__(self, structures=[], size=Vec3(1e-9), engine=None):
    super().__init__(engine)

    self.structures = structures if isinstance(structures, list) else [structures]
    self.default_sim_type = "bandstructure"
    self.size = size

  def _get_structures(self):
    return self.structures

  def _simulate_bandstructure(self, sess, nsweeps=20, fmax=2e15, krange=(0, 1), run_time=600e-15):
    frange = (F_EPSILON, fmax)

    sess.set_sim_region(size=self.size * Vec3(1.0, 1.3, 1.3))
    sess.set_sources(DipoleSource(frange=frange, pos=Vec3(0), axis=Vec3(0, 0, 1)))
    sess.set_sim_time(run_time)

    output = []
    for s in range(nsweeps):
      k = krange[0] + (s / (nsweeps - 1)) * (krange[1] - krange[0])
      sess.set_sim_region(boundaries={ "x": { "type": "bloch", "k": k } })

      sweep = sess.run(FindResonances(BBox(Vec3(0), self.size), frange, run_time / 2))
      output.append({
        "k": k,
        "freqs": list(map(lambda r : r[0], sweep))
      })
      logging.info("Sweep %d/%d completed" % (s + 1, nsweeps))
    
    return output

class Cavity1D(SimulationObject):
  def __init__(self, unit_cells=[], additional_structures=[], engine=None):
    self.additional_structures = additional_structures if isinstance(additional_structures, list) else [additional_structures]
    self.unit_cells = unit_cells if isinstance(unit_cells, list) else [unit_cells]

    super().__init__(engine)
    self.default_sim_type = "resonance"

  def _get_structures(self):
    if len(self.unit_cells) < 1:
      return []

    structs = []
    offset = Vec3(0)
    for c in self.unit_cells:
      for s in c.structures:
        copy_s = s.copy()
        copy_s.pos = Vec3(offset)
        copy_s.pos.x += c.size.x / 2
        structs.append(copy_s)

      offset.x += c.size.x

    for s in structs:
      s.pos -= offset / 2

    for s in self.additional_structures:
      copy_s = s.copy()
      s.material.mesh_order = 3
      structs.append(copy_s)

    return structs

  def _getsize(self):
    x = 0
    y = 0
    z = 0

    for c in self.unit_cells:
      x += c.size.x
      y = max(y, c.size.y)
      z = max(z, c.size.z)

    return Vec3(x, y, z)

  def _simulate_resonance(self, sess, target_freq=1e15, fspan=2e14, run_time=600e-15, analyze_time=590e-15):
    size = self._getsize()
    frange = (target_freq - fspan / 2, target_freq + fspan / 2)

    sess.set_sim_region(size=(size * 1.5))
    sess.set_sources(DipoleSource(frange=frange, pos=Vec3(0), axis=Vec3(0, 0, 1)))
    sess.set_sim_time(run_time)

    bbox = BBox(Vec3(0), size * 1.1)
    
    st = analyze_time
    et = st + 1 / target_freq

    sess.run({
      "res": FindResonances(bbox, frange, run_time / 2),
      "e": WaveEnergy(bbox, target_freq, st),
      "pxmin": SidePower(bbox, AXIS_X, -1, st, et),
      "pxmax": SidePower(bbox, AXIS_X, 1, st, et),
      "pymin": SidePower(bbox, AXIS_Y, -1, st, et),
      "pymax": SidePower(bbox, AXIS_Y, 1, st, et),
      "pzmin": SidePower(bbox, AXIS_Z, -1, st, et),
      "pzmax": SidePower(bbox, AXIS_Z, 1, st, et)
    }, analyze=False)

    res = sess.analyze("res")
    res = res[res[:,1].argsort()[0]]
    freq = res[0]

    e = sess.analyze("e", freq = freq)

    qxmin = 2*np.pi*freq*e/np.mean(sess.analyze("pxmin"))
    qxmax = 2*np.pi*freq*e/np.mean(sess.analyze("pxmax"))
    qymin = 2*np.pi*freq*e/np.mean(sess.analyze("pymin"))
    qymax = 2*np.pi*freq*e/np.mean(sess.analyze("pymax"))
    qzmin = 2*np.pi*freq*e/np.mean(sess.analyze("pzmin"))
    qzmax = 2*np.pi*freq*e/np.mean(sess.analyze("pzmax"))

    return {
      "freq": freq,
      "qxmin": qxmin,
      "qxmax": qxmax,
      "qymin": qymin,
      "qymax": qymax,
      "qzmin": qzmin,
      "qzmax": qzmax,
    }
