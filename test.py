from wvgsolver.simulation.objects import UnitCell, Cavity1D
from wvgsolver.geometry.structures import BoxStructure, CylinderStructure
from wvgsolver.utils.linalg import Vec3, BBox
from wvgsolver.geometry.materials import DielectricMaterial
import logging
import sys
import matplotlib.pyplot as plt
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

ncells = 6
w = 500e-9
h = 500e-9
r = 75e-9
nbeam = 2.4

cells = []
ta = 0
for i in range(ncells):
  a = 270e-9 - 30e-9 * (1 - ((i - ((ncells - 1) / 2)) / ((ncells - 1) / 2))**2)
  cells.append(UnitCell(
    structures=CylinderStructure(Vec3(0), h, r, DielectricMaterial(1.0)),
    size=Vec3(a, w, h)
  ))
  ta += a

cavity = Cavity1D(
  cells,
  additional_structures=BoxStructure(Vec3(0), Vec3(ta, w, h), DielectricMaterial(nbeam))
)
#r = cavity.simulate("resonance", target_freq=4.3e14, fspan=1e14)
#print(r)

cells[0].structures.append(
  BoxStructure(Vec3(0), Vec3(270e-9, w, h), DielectricMaterial(nbeam))
)
res = cells[0].simulate("bandstructure", fmax=1e15, nsweeps=60)

for r in res:
  plt.scatter([r["k"]] * len(r["freqs"]), r["freqs"], c="b")

plt.show()
