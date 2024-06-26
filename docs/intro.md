# Welcome to pyISFenix docs

**pyISFenix**: PYthon framework for the Incompressible Schrodinger Flow using FENIcsX

**Authors**: Stefano Riva, Carolina Introini, Antonio Cammi

The simulation of fluid dynamics is typically carried out either under the hypothesis of continuum matter (state-of-the-art Computational Fluid Dynamics) or adopting a particle perspective in Lattice Boltzmann methods; however, back in the early 1920s, while quantum mechanics theory was being developed, Erwin Madelung published a paper {cite:p}`Madelung1926` in which an analogy between the Schrödinger equation and the Euler equations, describing the motion of inviscid flows. Thus, a transfomation exist linking wave functions and flow fields.

This approach has been recently re-discovered by [Albert Chern](https://cseweb.ucsd.edu/~alchern/) in {cite:p}`Chern2017, Chern2016`, who published a numerical implementation of the Madelung transformation, enabling the simulation of incompressible inviscid fluids using an correspondent wave function. This novel methodology is known as Incompressible Schrödinger Flow (ISF): the original implementation has been later extended by the authors to a Finite Element solution by the authots in {cite:p}`Riva_ISF`; in particular, the [FEniCSx](https://fenicsproject.org/) library in Python is used to solve the steps of the ISF algorthm.

This documentation includes a brief introduction to the algorithm and the basis theory behind, the API documentation of the solvers and some case studies reported in the published papers.

This work has been carried out at the [Nuclear Reactors Group](https://www.nuclearenergy.polimi.it) at [Politecnico di Milano](https://polimi.it), under the supervision of Prof. Antonio Cammi.

---

## How to cite

If you are going to use this implementation of *ISF* in your research work, please cite the following articles

- Riva, S., Introini C., Cammi A., ‘A Finite Element Implementation of the Incompressible Schrödinger Flow Method’. Physics of Fluids, 2023.

---

## Installation notes

The package has been tested on MacOS and Ubuntu 18.04/20.04 machines with **Python3.10**.

### Dependencies
The *ISF* solvers requires the following dependencies:

```python
import numpy
import h5py

import pyvista
import gmsh
import dolfinx
```

Be sure to install *gmsh* and *gmsh-api* before *dolfinx=0.7.2* (the package has been tested with real mode of the PETSc library). The instructions to install *dolfinx* are available at [https://github.com/FEniCS/dolfinx#binary](https://github.com/FEniCS/dolfinx#binary).

### Set up a conda environment for *ISF*

At first create a new conda environment
```bash
conda create --name <env_name>
```
If not already done, add conda-forge to the channels
```bash
conda config --add channels conda-forge
```
After having activate it, install
```bash
conda install python=3.10
```
This provides also *pip* which is necessary to install *gmsh* as
```bash
python -m pip install gmsh gmsh-api
```
Now, install *dolfinx* (the current version works with version 0.7.2, an update to version 0.8.0 is under development)
```bash
conda install fenics-dolfinx=0.7.2 petsc=*=complex* mpich pyvista
```
Add the following packages
```bash
conda install meshio tqdm
```
Downgrade the following
```bash
python -m pip install setuptools==62.0.0
conda install numpy=1.23.5
```
Once this is completed, it may be necessary to re-install *gmsh*
```bash
python -m pip install gmsh gmsh-api
```

### How to use the solvers?

Once all the dependencies have been installed, the classes in `src` folder can be imported and used to solve ISF.
