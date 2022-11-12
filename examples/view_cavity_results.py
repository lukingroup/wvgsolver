"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial
from wvgsolver.engine import LumericalEngine
import os

#  Initialize Lumerical File Locations
FDTDLoc = '/n/home08/eknall/sw_ENK/lumerical-2021-R2-2717-7bf43e7149_seas'
FDTDexeLoc = os.path.join(FDTDLoc,'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc,'bin/fdtd-engine-ompi-lcl')

# Unit cells are cubes
cell_size = 0.3e-6
# Radius of the air holes in the cells
hole_radius = 0.075e-6
# The length of the cavity beam
beam_length = 10e-6
# The target resonance frequency, in Hz
target_frequency = 400e12

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=1, hide=True, lumerical_path=FDTDLoc, working_path="./fsps")

cell_box = BoxStructure(Vec3(0), Vec3(cell_size), DielectricMaterial(2, order=2, color="red"))
cell_hole = CylinderStructure(Vec3(0), cell_size, hole_radius, DielectricMaterial(1, order=1, color="blue"))

mirror_cells = [UnitCell(structures=[ cell_box, cell_hole ], size=Vec3(cell_size), engine=engine)] * 4
cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(cell_size), engine=engine)] * 3

cavity = Cavity1D(
  unit_cells=mirror_cells + cavity_cells + mirror_cells,
  structures=[ BoxStructure(Vec3(0), Vec3(beam_length, cell_size, cell_size), DielectricMaterial(2, order=2, color="red")) ],
  engine=engine
)

# By setting the save path here, the cavity will save itself after each simulation to this file
# cavity.save("cavity.obj")

r3 = cavity.simulate("guidedness", target_freq=target_frequency)
r1 = cavity.simulate("resonance", target_freq=target_frequency)

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
  r1["freq"], r1["vmode"],
  1/(1/r1["qxmin"] + 1/r1["qxmax"]),
  1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

r2 = cavity.simulate("quasipotential", target_freq=target_frequency)

# Plot the quasipotential
r2.show()


print("Guidedness: %s" % r3)
