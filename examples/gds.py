from wvgsolver import Cavity1D
from wvgsolver.parse import DielectricExtrusionFaceGDSParser
from wvgsolver.engine import LumericalEngine

engine = LumericalEngine()
cavity = Cavity1D(load_path="cavity.obj", engine=engine)

parsed = DielectricExtrusionFaceGDSParser(cavity)
parsed.show()
