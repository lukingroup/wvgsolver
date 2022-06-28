"""
This example loads a cavity from a file, and inspects its results. The same
thing can be done for any simulation object
"""

from wvgsolver import Cavity1D
from wvgsolver.engine import LumericalEngine
import datetime

engine = LumericalEngine()

cavity = Cavity1D(load_path="cavity.obj", engine=engine)

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
