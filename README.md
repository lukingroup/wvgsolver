# Waveguide and Unit Cell Solver Package

## Installation

```
pip install .
```

From the root repository directory

## Documentation

This package allows for running simulations on 1 dimension nanophotonic cavities
and their constituent unit cells. For these two objects, we expose two classes:
`Cavity1D` and `UnitCell`, which derive from the more generic `SimulationObject`, and
provide functionality to run simulations and store and visualize results.
3D geometry for these objects is represented by Structure and Material classes.
Configuring your simulation environment can be done using Engine classes, and
cavities and unit cells can be parsed into GDS files using Parser classes.

### Creating an object

When creating a new simulation object, you pass the parameters that describe its 3D structure. To create a
unit cell, simply pass a list of structures in its constructor. For example, creating a unit
cell consisting of a dielectric box with an air hole can be done as follows:

```
cell = UnitCell(
  structures=[
    CylinderStructure(Vec3(0), 1e-6, 1e-6, DielectricMaterial(1, order=1)),
    BoxStructure(Vec3(0), Vec3(1e-6, 1e-6, 1e-6), DielectricMaterial(2, order=2))
  ],
  size=Vec3(1e-6, 1e-6, 1e-6)
)
```

There are three available types of structures right now: An extruded polygon, a box, and
an elliptical cylinder. Each structure also has an associated material, such as a dielectric.

To create a cavity, pass a list of unit cells plus any additional structures (such as a beam). For example,

```
cavity = Cavity1D([cell], structures=[
  BoxStructure(Vec3(0), Vec3(10e-6, 1e-6, 1e-6), DielectricMaterial(2, order=2))
])
```

See `wvgsolver/simulation/objects.py` for a description of the constructor parameters of each
simulation object class, and `wvgsolver/geometry/structures.py` for documentation on how to define
structures.

### Simulating objects

Once you've created an object, you can call its various simulation methods to get some data! This is done
using the `simulate` method, which takes as its first parameter the type of simulation, as well
as additional keyword parameters specific to that simulation type.

With a `Cavity1D` object, you can simulate its resonance frequency:

```
cavity.simulate("resonance", target_freq=407e12)
```

Or the quasipotential of its cells:

```
cavity.simulate("quasipotential", target_freq=407e12)
```

Or its guidedness:

```
cavity.simulate("guidedness", target_freq=407e12)
```

With a `UnitCell` object, you can simulate its bandstructure:

```
cell.simulate("bandstructure", ks=(0, 0.5, 50), freqs=(200e12, 600e12, 10000))
```

Or its bandgap:

```
cell.simulate("bandgap", freqs=(200e12, 600e12, 10000))
```

Each simulate call will return the results from that simulation.
The full list of parameters, and the return type, for each simulation type
is described in `wvgsolver/simulation/objects.py`.

### Getting results

In addition to returning its results, every simulation call caches its results in its parent
object's internal data. You can then access this, optionally filtering for results
from simulations that used particular parameters, using the `get_results()` method. For
example, to get all results from "resonance" simulations on a cavity, call

```
res = cavity.get_results("resonance")
```

See `wvgsolver/simulation/base.py` for a description of the parameters of `get_results()`.

### Saving and loading objects

Simulation objects can be saved to files by calling `save()`. For example, you can save
a cavity with

```
cavity.save("cavity.obj")
```

After saving an object for the first time, after every subsequent simulation the object
will push its new list of simulation results to that file, unless you pass `save=False`
to `simulate`.

To load an object that you previously saved to a file, use the `load_path` parameter:

```
cavity = Cavity1D(load_path="cavity.obj")
```

### Customizing the Engine

By default, all simulations use Lumerical FDTD. You can customize some parameters
of your Lumerical simulation environment by creating an instance of LumericalEngine,
and passing it to the objects you intend to simulate. For example, to create an engine
that uses an especially coarse mesh, and stores its files in a particular directory:

```
engine = LumericalEngine(working_path="fsp_files/", mesh_accuracy=1)
cavity = Cavity1D([cell], engine=engine)
```

See `wvgsolver/engine/engines.py` for more information on LumericalEngine options.
