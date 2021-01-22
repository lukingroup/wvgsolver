from .base import Structure

class BoxStructure(Structure):
  def __init__(self, pos, size, material):
    super().__init__(pos)
    self.size = size
    self.material = material

  def _add_lumerical(self, session):
    self.material.add(session)
    session.fdtd.addrect(name=self.name, x=self.pos.x, y=self.pos.y, z=self.pos.z,
      x_span=self.size.x, y_span=self.size.y, z_span=self.size.z, material=self.material.name)

class CylinderStructure(Structure):
  def __init__(self, pos, height, radius, material):
    super().__init__(pos)
    self.height = height
    self.radius = radius
    self.material = material

  def _add_lumerical(self, session):
    self.material.add(session)
    session.fdtd.addcircle(name=self.name, x=self.pos.x, y=self.pos.y, z=self.pos.z,
      z_span=self.height, radius=self.radius, material=self.material.name)
