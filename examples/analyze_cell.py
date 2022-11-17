"""
This example loads a cavity from a file, and inspects its results. The same
thing can be done for any simulation object
"""

from wvgsolver import UnitCell
from wvgsolver.engine import LumericalEngine
import datetime
import sys
import os

#  Initialize Lumerical File Locations
FDTDLoc = '/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas'

engine = LumericalEngine(lumerical_path=FDTDLoc)

cell_name = sys.argv[1]
cell = UnitCell(load_path=cell_name, engine=engine)

# Get all resonance simulation results whose simulations started after 4/20/2021
res = cell.get_results("bandstructure", started_after=datetime.datetime(2021, 4, 20))

# Filter any failed results out
res = list(filter(lambda r : (r["status"] == "succeeded"), res))

if len(res):
  # Take the latest result
  r = res[-1]

  # Let's open the Lumerical FSP file used for this simulations
  print(r)
  print(r["sess_res"])
  r["sess_res"].show()
  #r["res"]["xyprofile"].show()
  #r["res"]["yzprofile"].show()


  #r["res"]["xyprofile"].save("temp_xy.png")
