from distutils.core import setup

setup(
  name = "wvgsolver",
  packages = ["wvgsolver", "wvgsolver.simulation", "wvgsolver.parse",
    "wvgsolver.geometry", "wvgsolver.engine", "wvgsolver.utils"],
  py_modules=["wvgsolver.analysis.base", "wvgsolver.analysis.procedures"],
  version = "0.1",
  author = "Vassilios Kaxiras",
  author_email = "vassilioskaxiras@gmail.com",
  description = "1D Cavity Waveguide and Unit Cell simulations",
)
