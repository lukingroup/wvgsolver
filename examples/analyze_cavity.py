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
FDTDLoc = '/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas'

engine = LumericalEngine(lumerical_path=FDTDLoc)

cavity_name = sys.argv[1]
cavity = Cavity1D(load_path=cavity_name, engine=engine)


r2 = cavity.simulate("quasipotential", target_freq=406.7e12)

# Plot the quasipotential
r2.show()

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
  r["res"]["xyprofile"].show()
  print(r["res"]["xyprofile"])
  print(type(r["res"]["xyprofile"]))
  r["res"]["xyprofile"].save("temp_xyprofile.png",title="title check")
  r["res"]["yzprofile"].show()


  r["res"]["xyprofile"].save("temp_xy.png")
