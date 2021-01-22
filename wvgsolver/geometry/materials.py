from .base import Material

class DielectricMaterial(Material):
  def __init__(self, nindex):
    super().__init__()
    self.nindex = nindex

  def _add_lumerical(self, sess):
    if sess.fdtd.materialexists(self.name):
      sess.fdtd.deletematerial(self.name)
    mat = sess.fdtd.addmaterial("Dielectric")
    sess.fdtd.setmaterial(mat, "name", self.name)
    sess.fdtd.setmaterial(self.name, "Refractive Index", self.nindex)
