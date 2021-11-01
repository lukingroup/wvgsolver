from .base import Material

class DielectricMaterial(Material):
  """Represents a Dielectric material for a structure"""
  def __init__(self, nindex, order=2, color="red"):
    """
    Parameters
    ----------
    nindex : float
      The index of refraction of the material
    order : int
      The priority order of the structure. Lower means higher priority, down to 0. In regions
      where structures overlap, the structure with the lower order override the one with the
      higher order
    color : Matplotlib color
      The color of the material. This is used by some engines when displaying the structure.
      Can be a string such as "red", an array of RGB or RGBA values, and other values.
    """
    super().__init__(order, color)
    self.nindex = nindex

  def __repr__(self):
    return "DielectricMaterial(%s):%s" % (self.nindex, super().__repr__())

  def _add_lumerical(self, sess):
    if sess.fdtd.materialexists(self.name):
      sess.fdtd.deletematerial(self.name)
    mat = sess.fdtd.addmaterial("Dielectric")
    sess.fdtd.setmaterial(mat, "name", self.name)
    sess.fdtd.setmaterial(self.name, "Refractive Index", self.nindex)
    sess.fdtd.setmaterial(self.name, "Mesh Order", self.order)

    sess.fdtd.setmaterial(self.name, "color", self.get_color_rgba())
