import math

# TODO: Double check this function
def axis_to_spherical(uaxis):
  """Convert an axis to spherical coordinates"""
  axis = uaxis.normalized()
  theta = math.acos(axis.z) * 180 / math.pi
  phi = math.copysign(90, axis.y)
  if axis.x != 0.0:
    phi = math.atan(axis.y/axis.x) * 180 / math.pi

  return theta, phi

class Vec3:
  """Simple class for 3-vector linear algebra"""
  def __init__(self, x, y=None, z=None):
    if isinstance(x, Vec3):
      self.x = x.x
      self.y = x.y
      self.z = x.z
    else:
      self.x = x
      if y is None:
        self.y = x
      else:
        self.y = y
      if z is None:
        self.z = x
      else:
        self.z = z

  def __repr__(self):
    return "Vec3(%.6e,%.6e,%.6e)" % (self.x, self.y, self.z)

  def __hash__(self):
    return hash(self.__repr__())

  def __eq__(self, other):
    if isinstance(other, Vec3):
      return self.x == other.x and self.y == other.y and self.z == other.z
    return False

  def __add__(self, other):
    if isinstance(other, Vec3):
      return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
    
    return Vec3(self.x + other, self.y + other, self.z + other)
  
  def __iadd__(self, other):
    if isinstance(other, Vec3):
      self.x += other.x
      self.y += other.y
      self.z += other.z
    else:
      self.x += other
      self.y += other
      self.z += other

    return self
  
  def __sub__(self, other):
    if isinstance(other, Vec3):
      return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    return Vec3(self.x - other, self.y - other, self.z - other)
  
  def __isub__(self, other):
    if isinstance(other, Vec3):
      self.x -= other.x
      self.y -= other.y
      self.z -= other.z
    else:
      self.x -= other
      self.y -= other
      self.z -= other

    return self
  
  def __floordiv__(self, other):
    if isinstance(other, Vec3):
      return Vec3(self.x // other.x, self.y // other.y, self.z // other.z)
    
    return Vec3(self.x // other, self.y // other, self.z // other)
  
  def __ifloordiv__(self, other):
    if isinstance(other, Vec3):
      self.x //= other.x
      self.y //= other.y
      self.z //= other.z
    else:
      self.x //= other
      self.y //= other
      self.z //= other

    return self
  
  def __truediv__(self, other):
    if isinstance(other, Vec3):
      return Vec3(self.x / other.x, self.y / other.y, self.z / other.z)
    
    return Vec3(self.x / other, self.y / other, self.z / other)
  
  def __itruediv__(self, other):
    if isinstance(other, Vec3):
      self.x /= other.x
      self.y /= other.y
      self.z /= other.z
    else:
      self.x /= other
      self.y /= other
      self.z /= other

    return self
  
  def __mul__(self, other):
    if isinstance(other, Vec3):
      return Vec3(self.x * other.x, self.y * other.y, self.z * other.z)
    
    return Vec3(self.x * other, self.y * other, self.z * other)
  
  def __imul__(self, other):
    if isinstance(other, Vec3):
      self.x *= other.x
      self.y *= other.y
      self.z *= other.z
    else:
      self.x *= other
      self.y *= other
      self.z *= other

    return self

  def length(self):
    return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

  def normalized(self):
    r = self.length()
    return Vec3(self.x / r, self.y / r, self.z / r)
  
  def tolist(self):
    return [self.x, self.y, self.z]

class BBox:
  """Simple class storing a bounding box"""
  def __init__(self, pos, size):
    """
    Parameters
    ----------
    pos : Vec3
      The position of the centroid of the box
    size : Vec3
      The size of each size of the box
    """
    self.pos = pos
    self.size = size
