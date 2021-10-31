from wvgsolver import Cavity1D
from wvgsolver.parse import DielectricExtrusionFaceGDSParser

cavity = Cavity1D(load_path="cavity.obj")

parsed = DielectricExtrusionFaceGDSParser(cavity, invert=True)
parsed.show()
