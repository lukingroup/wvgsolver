"""
This example creates a unit cell, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial
from wvgsolver.engine import LumericalEngine
import numpy as np


#  Initialize Lumerical File Location
FDTDLoc = '/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas'

# Unit cells are rectangular prisms
a = 0.215e-6
beam_w = 0.45*1e-6
beam_h = 0.2*1e-6
# Dimensions of the air hole in the unit cell
hx = 0.5228*a
hy = 0.7265*a


# The target resonance frequency, in Hz
target_frequency = 406.7e12

# Use level 3 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=3, hide=False, lumerical_path=FDTDLoc)
cell_size = Vec3(a,beam_w,beam_h)
cell_box = BoxStructure(Vec3(0), cell_size, DielectricMaterial(2.4028, order=2, color="blue"))

# offset the hole to respect the way we define the relevant lattice constant
cell_hole = CylinderStructure(Vec3(-a/2,0,0), beam_h, hx/2, DielectricMaterial(1, order=1, color="red"), radius2=hy/2)
cell = UnitCell(structures = [cell_box, cell_hole], size=cell_size, engine=engine)





# By setting the save path here, the cavity will save itself after each simulation to this file
cell.save("cell.obj")

r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 12), freqs=(0.25e15, 0.7e15, 100000), 
                    dipole_region=Vec3(0.8, 0, 0), window_pos = 0)

# Plot the bandstructure
r1.show()


#r2 = cell.simulate("bandgap")

#print("Bandgap: %f Hz" % (r2[1] - r2[0]))
