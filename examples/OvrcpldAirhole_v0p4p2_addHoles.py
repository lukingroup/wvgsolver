"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial
from wvgsolver.engine import LumericalEngine

import numpy as np

# Unit cells are triangular prisms
beam_width = 0.480e-6
apex_half_angle = 50*np.pi/180
lattice_constant = 0.270e-6
# Radius of the air holes in the cells
hole_radius = 0.5*0.120e-6
# The length of the cavity beam
beam_length = 10e-6
# The target resonance frequency, in Hz
target_frequency = 400e12

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=4, hide=False)

cell_box = TriStructure(Vec3(0), Vec3(beam_width, apex_half_angle, lattice_constant), DielectricMaterial(2.43, order=2))
cell_hole = CylinderStructure(Vec3(0), lattice_constant, hole_radius, DielectricMaterial(1, order=1))

mirror_cells = [UnitCell(structures=[ cell_box, cell_hole ], size=Vec3(lattice_constant), engine=engine)] * 4
cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(lattice_constant), engine=engine)] * 3

cavity = Cavity1D(
  unit_cells=mirror_cells + cavity_cells + mirror_cells,
  structures=[ BoxStructure(Vec3(0), Vec3(beam_length, beam_width, beam_width), DielectricMaterial(2, order=2)) ],
  engine=engine
)

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("v0p4p2p2.obj")

r1 = cavity.simulate("resonance", target_freq=target_frequency)

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
  r1["freq"], r1["vmode"],
  1/(1/r1["qxmin"] + 1/r1["qxmax"]),
  1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

# r2 = cavity.simulate("quasipotential", target_freq=target_frequency)

# # Plot the quasipotential
# r2.show()

# r3 = cavity.simulate("guidedness", target_freq=target_frequency)

# print("Guidedness: %f" % r3)

