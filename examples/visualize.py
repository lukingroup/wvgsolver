from wvgsolver import Cavity1D
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D

cavity = Cavity1D(load_path="cavity.obj")

parsed = ObjectParser3D(cavity)
parsed.show()

parsed = DielectricExtrusionFaceGDSParser(cavity, invert=True)
parsed.show()
