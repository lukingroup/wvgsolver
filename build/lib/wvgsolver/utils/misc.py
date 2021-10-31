import uuid

def hasmethod(obj, name):
  return hasattr(obj, name) and callable(getattr(obj, name))

def randstring():
  return str(uuid.uuid4())
