from .base import Engine, Session
try:
  import meep as mp
except:
  pass
from ..utils.meep import U_A, U_K, U_F, U_T
from ..utils.constants import AXIS_X, AXIS_Y, AXIS_Z
import numpy as np
from scipy.sparse.linalg import eigs
from scipy.sparse import diags, eye, bmat
import matplotlib.pyplot as plt
import trimesh

class EffIndex1DSession(Session):
  def __init__(self, engine):
    super().__init__(engine)
    self.structures = []
    self.sources = []
    self.cell_size = mp.Vector3(1)
    self.boundary_layers = [mp.PML(self.engine.pml_thickness/U_A)]
    self.sim_time = 1
    self.sim = None
    self.analyze_regions = []
    self.symmetric_axes = []

  def _set_structures(self, structures):
    self.structures = structures

  def _set_sources(self, sources):
    self.sources = [ s.add(self) for s in sources ]
  
  def _set_mesh_regions(self, regions=[]):
    pass

  def add_analyze_region(self, region):
    self.analyze_regions.append(region)
  
  def remove_analyze_region(self, region):
    self.analyze_regions.remove(region)

  def close(self):
    pass

  def set_sim_region(self, pos=None, size=None, boundaries={}):
    if size is not None:
      self.cell_size = mp.Vector3(size.x/U_A, size.y/U_A, size.z/U_A)

    axis_map = {
      "x": AXIS_X,
      "y": AXIS_Y,
      "z": AXIS_Z
    }

    for a in ["x", "y", "z"]:
      for s in ["min", "max"]:
        for t in ["antisymmetric", "symmetric"]:
          if (a + s) in boundaries and boundaries[a + s] == t:
            axis = axis_map[a]
            if axis not in self.symmetric_axes:
              self.symmetric_axes.append(axis)

  def set_sim_time(self, t):
    self.sim_time = t/U_T

  def _compute_mode(self, f, midx, x, y, eps):
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    X, Y = np.meshgrid(x, y)

    k0 = 2*np.pi*f
    eps = eps.flatten()

    Uy = diags([-1, 1], [0, 1], (eps.shape[0], eps.shape[0]))/dy
    Ux = diags([-1, 1], [0, y.shape[0]], (eps.shape[0], eps.shape[0]))/dx
    Vy = -np.transpose(Uy)
    Vx = -np.transpose(Ux)

    vepsrx = (np.concatenate(([1], eps[:-1])) + eps)/2
    vepsry = (np.concatenate((np.ones(y.shape[0]), eps[:-y.shape[0]])) + eps)/2
    vepsrz = (np.concatenate((np.ones(y.shape[0]), eps[:-y.shape[0]])) + \
             np.concatenate((np.ones(y.shape[0] + 1), eps[:-(y.shape[0] + 1)])) + \
             np.concatenate(([1], eps[:-1])) + eps)/4

    epsrx = diags(vepsrx)
    epsry = diags(vepsry)
    epsrz = diags(vepsrz)
    iepsrx = diags(1/vepsrx)
    iepsry = diags(1/vepsry)
    iepsrz = diags(1/vepsrz)

    I = eye(eps.shape[0]) 
    
    Pxx = -(1/k0**2)*Ux @ iepsrz @ Vy @ Vx @ Uy + (k0**2*I + Ux @ iepsrz @ Vx) @ (epsrx + (1/k0**2)*Vy @ Uy)
    Pyy = -(1/k0**2)*Uy @ iepsrz @ Vx @ Vy @ Ux + (k0**2*I + Uy @ iepsrz @ Vy) @ (epsry + (1/k0**2)*Vx @ Ux)
    Pxy = Ux @ iepsrz @ Vy @ (epsry + (1/k0**2)*Vx @ Ux) - (1/k0**2)*(k0**2*I + Ux @ iepsrz @ Vx) @ Vy @ Ux
    Pyx = Uy @ iepsrz @ Vx @ (epsrx + (1/k0**2)*Vy @ Uy) - (1/k0**2)*(k0**2*I + Uy @ iepsrz @ Vy) @ Vx @ Uy

    P = bmat([[Pxx, Pxy], [Pyx, Pyy]])

    b2, s = eigs(P, k=(midx+1), which="LR")

    neff = np.sqrt(b2[-1]) / k0
    s = np.abs(s[:,-1].reshape((2, x.shape[0], y.shape[0])))**2
    s = np.sum(s, axis=0) / np.sum(s)
    return neff, s
    
  def _compute_epsilons(self):
    nmesh = self.cell_size / (self.engine.resolution/U_A)

    npointsx = int(nmesh.x) if AXIS_X not in self.symmetric_axes else int(nmesh.x / 2)
    npointsy = int(nmesh.y) if AXIS_Y not in self.symmetric_axes else int(nmesh.y / 2)
    npointsz = int(nmesh.z) if AXIS_Z not in self.symmetric_axes else int(nmesh.z / 2)

    x = np.linspace(
      -self.cell_size.x/2 if AXIS_X not in self.symmetric_axes else 0,
      self.cell_size.x/2, npointsx
    )
    y = np.linspace(
      -self.cell_size.y/2 if AXIS_Y not in self.symmetric_axes else 0,
      self.cell_size.y/2, npointsy
    )
    z = np.linspace(
      -self.cell_size.z/2 if AXIS_Z not in self.symmetric_axes else 0,
      self.cell_size.z/2, npointsz
    )
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    npoints = npointsx*npointsy*npointsz
    eps = np.ones(npoints)
    points = np.stack((X, Y, Z), axis=-1).reshape((npoints, 3))

    for s in sorted(self.structures, key=lambda st: st.material.order, reverse=True):
      eps[s.contains(points)] = s.material.nindex**2

    eps = eps.reshape((npointsx, npointsy, npointsz))
    grid_axes = [x, y, z]
    for i, a in enumerate([AXIS_X, AXIS_Y, AXIS_Z]):
      if a in self.symmetric_axes:
        indexer = [np.s_[:], np.s_[:], np.s_[:]]
        indexer[i] = np.s_[:-1]
        eps = np.concatenate((np.flip(eps, axis=i)[tuple(indexer)], eps), axis=i)
        grid_axes[i] = np.concatenate((-np.flip(grid_axes[i])[:-1], grid_axes[i]))

    eps_r = eps[np.abs(grid_axes[0] - self.engine.reference_point/U_A).argmin(),:,:]
    neff, mode = self._compute_mode(self.engine.mode_f/U_F, self.engine.mode_index, grid_axes[1], grid_axes[2], eps_r)
    self.engine.last_mode = neff, np.flip(mode.T, axis=0)
  
    return np.real(neff**2 + np.sum((eps - np.expand_dims(eps_r, 0))*np.expand_dims(mode, 0), axis=(1, 2)))

  def _prerun(self):
    epsilons = self._compute_epsilons()

    geometry = [
      mp.Block(
        mp.Vector3(mp.inf, mp.inf, self.cell_size.x),
        center=mp.Vector3(),
        material=mp.MaterialGrid(
          (1, 1, epsilons.shape[0]),
          mp.Medium(epsilon=np.min(epsilons)),
          mp.Medium(epsilon=np.max(epsilons)),
          (epsilons - np.min(epsilons)) / (np.max(epsilons) - np.min(epsilons))
        )
      )
    ]

    self.sim = mp.Simulation(
      cell_size=mp.Vector3(z=self.cell_size.x),
      dimensions=1,
      boundary_layers=self.boundary_layers,
      geometry=geometry,
      sources=self.sources,
      resolution=(1/self.engine.resolution)*U_A
    )

  def _runsim(self):
    regions = sorted(self.analyze_regions, key=lambda r: r["time"])
    last_t = 0
    for r in regions:
      t = r["time"]
      self.sim.run(until=(t - last_t))
      obj = None
      if r["type"] == "flux":
        obj = self.sim.add_flux(*r["args"], **r["kwargs"])
      elif r["type"] == "energy":
        obj = self.sim.add_energy(*r["args"], **r["kwargs"])
      last_t = t
      r["callback"](obj)

    self.sim.run(until=(self.sim_time - last_t))

    sim_region = mp.Vector3(z=self.cell_size.x)
    def debug(sim):
      e_data = sim.get_array(center=mp.Vector3(), size=sim_region, component=mp.Ex)
      eps_data = sim.get_array(center=mp.Vector3(), size=sim_region, component=mp.Dielectric)
      plt.clf()
      plt.plot(np.real(eps_data), label="Re[epsilon]")
      plt.plot(100*np.abs(e_data), label="Ex")
      plt.legend()
      plt.pause(0.1)

#    self.sim.run(mp.at_every(1, debug), until=0.1*self.sim_time)

"""
Defines a 1d simulation engine that reduces a 3d structure to 1d along the x axis using effective index methods
"""
class EffIndex1DEngine(Engine):
  def __init__(self, resolution=10e-9, pml_thickness=1e-6, reference_point=0, mode_index=0, mode_f=407e12):
    self.name = "eff1d"
    self.resolution = resolution
    self.pml_thickness = pml_thickness
    self.reference_point = reference_point
    self.mode_index = mode_index
    self.mode_f = mode_f
    self.last_mode = None

  def new_session(self):
    return EffIndex1DSession(self)
