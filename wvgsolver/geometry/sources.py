from .base import Source
from ..utils.constants import C_LIGHT, AXIS_Z
from ..utils.linalg import Vec3
import math

class DipoleSource(Source):
  def __init__(self, lrange=None, frange=None, pos=Vec3(0.0), axis=Vec3(0.0, 0.0, 1.0)):
    super().__init__(self)
    if lrange is None and frange is None:
      raise TypeError("One of lrange and frange must be set")

    if lrange is not None:
      self.lambda_start = lrange[0]
      self.lambda_end = lrange[1]
    elif frange is not None:
      self.lambda_start = C_LIGHT / frange[1]
      self.lambda_end = C_LIGHT / frange[0]

    self.axis = axis
    self.pos = pos

  def _add_lumerical(self, session):
    axis = self.axis.normalized()
    theta = math.acos(axis.z) * 180 / math.pi
    phi = math.copysign(90, axis.y)
    if axis.x != 0.0:
      phi = math.atan(axis.y/axis.x) * 180 / math.pi

    session.fdtd.adddipole(x=self.pos.x, y=self.pos.y, z=self.pos.z,
      wavelength_start=self.lambda_start, wavelength_stop=self.lambda_end,
      theta=theta, phi=phi, name=self.name)

