from wvgsolver import Cavity1D
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D

cavity = Cavity1D(load_path="cavity.obj")

# Display the Lumerical file associated with a simulation
# cavity.get_results("guidedness")[0]["sess_res"].show()

parsed = ObjectParser3D(cavity)
# Display a 3D visualization of the cavity
parsed.show()

parsed = DielectricExtrusionFaceGDSParser(cavity, invert=True)
# Display a 2D slice of the cavity, similarly to what we might output to a GDS file
parsed.show()
