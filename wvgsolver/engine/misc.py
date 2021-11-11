class LumericalMeshRegion:
  """Defines a mesh override region in Lumerical FDTD"""
  def __init__(self, bbox, dx, dy=None, dz=None):
    """
    Parameters
    ----------
    bbox : BBox
      The bbox of the region
    dx : float
      The maximum mesh step, in meters. To specificy the same step in all directions, 
      only provide this value. If dy and dz are specified, this is only taken to be
      the step size in the x direction.
    dy : float or None
      The maximum mesh step, in meters, in the y direction.
    dz : float or None
      The maximum mesh step, in meters, in the z direction.
    """
    self.bbox = bbox
    self.dx = dx
    self.dy = dx
    self.dz = dx
    if dy is not None:
      self.dy = dy
    if dz is not None:
      self.dz = dz 

  def add(self, engine):
    engine.fdtd.addmesh(x=self.bbox.pos.x, y=self.bbox.pos.y, z=self.bbox.pos.z, \
      x_span=self.bbox.size.x, y_span=self.bbox.size.y, z_span=self.bbox.size.z, \
      dx=(self.dx*1e6), dy=(self.dy*1e6), dz=(self.dz*1e6))
