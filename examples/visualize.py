from wvgsolver import Cavity1D
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D

cavity = Cavity1D(load_path="cavity.obj")

cavity.get_results("guidedness")[0]["sess_res"].show()

parsed = ObjectParser3D(cavity)
parsed.show()

parsed = DielectricExtrusionFaceGDSParser(cavity, invert=True)
parsed.show()
