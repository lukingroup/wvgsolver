import math

class Vec3:
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

class BBox:
  def __init__(self, pos, size):
    self.pos = pos
    self.size = size
