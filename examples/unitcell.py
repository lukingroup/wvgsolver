"""
This example creates a unit cell, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial

# Unit cells are cubes
cell_size = 0.2e-6
# Radius of the air holes in the cells
hole_radius = 0.05e-6
# The target resonance frequency, in Hz
target_frequency = 400e12

# Use level 3 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=3, hide=False)

cell_box = BoxStructure(Vec3(0), Vec3(cell_size), DielectricMaterial(2, order=2)),
cell_hole = CylinderStructure(Vec3(0), cell_size, hole_radius, DielectricMaterial(1, order=1))

cell = UnitCell(structures=[ cell_box, cell_hole ])

# By setting the save path here, the cavity will save itself after each simulation to this file
cell.save("cell.obj")

r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 40))

# Plot the bandstructure
r1.show()

r2 = cell.simulate("bandgap")

print("Bandgap: %f Hz" % (r2[1] - r2[0]))