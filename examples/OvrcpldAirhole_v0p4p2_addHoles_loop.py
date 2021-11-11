"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D


import numpy as np
import os

base_dir = 'examples/v0p4p2/'
hole_radii_fns = ['holeStruct_4-4.txt',
                  'holeStruct_4-4.txt','holeStruct_3-3.txt',
                  'holeStruct_3-3.txt']
lat_const_fns = ['periodStruct_4-4_wrong.txt',
                'periodStruct_4-4.txt','periodStruct_3-3_wrong.txt',
                'periodStruct_3-3.txt']

beam_width = 0.482e-6
apex_half_angle = 50*np.pi/180
# lattice_constant = 0.270e-6
beam_height = (beam_width / 2) * np.tan(np.pi/2 - apex_half_angle)
# Radius of the air holes in the cells
# hole_radius = 0.5*0.120e-6
# The length of the cavity beam
beam_length = 15e-6
# The target resonance frequency, in Hz
target_frequency = 406.774e12

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=8, hide=False)

for h_fn,a_fn in zip(hole_radii_fns,lat_const_fns):

  hole_radii = np.loadtxt(os.path.join(os.path.curdir, base_dir+h_fn), dtype=float, usecols=(0,1), unpack=False)
  hole_radii /= 2
  lattice_constants = np.loadtxt(os.path.join(os.path.curdir, base_dir+a_fn), dtype=float, usecols=(0), unpack=False)
  # Unit cells are triangular prisms
  

  unit_cells = []
  for ((radius_x,radius_y), lattice_constant) in zip(hole_radii,lattice_constants):
    cell_box = TriStructure(Vec3(0), Vec3(beam_width, apex_half_angle, lattice_constant), 
                          DielectricMaterial(2.4028, order=2), rot_angles=(np.pi/2, np.pi/2, 0))
    
    cell_hole = CylinderStructure(Vec3(0), beam_height, radius_x, DielectricMaterial(1, order=1), radius2=radius_y)
    unit_cells += [UnitCell(structures=[ cell_box, cell_hole ], size=Vec3(lattice_constant,beam_width,beam_height), engine=engine)]

  cavity = Cavity1D(
    unit_cells=unit_cells,
    structures=[TriStructure(Vec3(0), Vec3(beam_width, apex_half_angle, beam_length), 
                          DielectricMaterial(2.4028, order=2), rot_angles=(np.pi/2, np.pi/2, 0))],
    engine=engine,
  )

  # parsed = ObjectParser3D(cavity)
  # parsed.show()

  # By setting the save path here, the cavity will save itself after each simulation to this file
  mirr_tag = a_fn[13:-4]
  cavity.save("v0p4p2p2_"+mirr_tag+".obj")

  r1 = cavity.simulate("resonance", target_freq=target_frequency)

  # Print the reults and plot the electric field profiles
  print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
    r1["freq"], r1["vmode"],
    1/(1/r1["qxmin"] + 1/r1["qxmax"]),
    1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
  ))

#   r1["xyprofile"].show()
#   r1["yzprofile"].show()

  cavity = Cavity1D(load_path="v0p4p2p2_"+mirr_tag+".obj",engine=engine)
  r1 = cavity.get_results("resonance")[0]
  print(r1['res']["xyprofile"].max_loc())
  print(r1['res']["yzprofile"].max_loc())
#   r1["sess_res"].show()

  print(r1)
r1["sess_res"].show()


  # print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
  #   r1["freq"], r1["vmode"],
  #   1/(1/r1["qxmin"] + 1/r1["qxmax"]),
  #   1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
  # ))

  # r2 = cavity.simulate("quasipotential", target_freq=target_frequency)

  # # Plot the quasipotential
  # r2.show()

  # r3 = cavity.simulate("guidedness", target_freq=target_frequency)

  # print("Guidedness: %f" % r3)

