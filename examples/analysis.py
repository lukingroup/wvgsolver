"""
This example loads a cavity from a file, and inspects its results. The same
thing can be done for any simulation object
"""

from wvgsolver import Cavity1D
from wvgsolver.engine import LumericalEngine
import datetime
import sys
import os

#  Initialize Lumerical File Locations
FDTDLoc = '/n/home08/eknall/sw_ENK/lumerical-2021-R2-2717-7bf43e7149_seas'
FDTDexeLoc = os.path.join(FDTDLoc,'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc,'bin/fdtd-engine-ompi-lcl')


engine = LumericalEngine(lumerical_path=FDTDLoc)

cavity_name = sys.argv[1]
cavity = Cavity1D(load_path=cavity_name, engine=engine)

# Get all resonance simulation results whose simulations started after 4/20/2021
res = cavity.get_results("resonance", started_after=datetime.datetime(2021, 4, 20))

# Filter any failed results out
res = list(filter(lambda r : (r["status"] == "succeeded"), res))

if len(res):
  # Take the latest result
  r = res[-1]

  # Print the resonance frequency and the keyword arguments passed to `simulate`
  # for this simulation
  print("Resonance freq: %f, Simulation kwargs: %s" % (
    r["res"]["freq"],
    r["kwargs"]
  ))

  # Let's open the Lumerical FSP file used for this simulation
  r["sess_res"].show()
  r["xyprofile"].show()
  r["yzprofile"].show()
