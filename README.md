# MP-pyISFenix

[![Original Paper [PoF]](https://img.shields.io/badge/Original%20Paper%20%5BPoF%5D-10.1063/5.0180356-gray?labelColor=blue&style=flat&link=https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)](https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)

**PYthon framework for the Incompressible Schrodinger Flow using FENIcsX**

**Authors:** Stefano Riva, Carolina Introini, Antonio Cammi

This repository collects some codes implemented to simulate inviscid fluids by means of the ISF technique, firstly developed by [Albert Chern](https://cseweb.ucsd.edu/~alchern/), adopting Finite Element Methods exploiting the [FEniCSx Library (v. 0.7.2)](https://fenicsproject.org/) for Python, extension to v0.8.0 is under development.

A **webpage** has been created at the following link [https://ermete-lab.github.io/MP-pyISFenix/intro.html](https://ermete-lab.github.io/MP-pyISFenix/intro.html), including the API documentation of the solvers and some tutorials resembling the case studies of the published papers.

--------------------------------

The ISF algorithm is based on the Madelung transform, an analogy between the wave functions governed by the Schrodinger equation and the velocity field for inviscid fluids.

Here's some useful background references:
- E. Madelung, Eine anschauliche deutung der gleichung von schr ̈odinger, Naturwissenschaften 14 (1926) 1004. doi:10.1007/BF01504657.
- A. Chern. Fluid Dynamics with Incompressible Schroedinger Flow. PhD thesis, 2017
- A. Chern, F. Knöppel, U. Pinkall, P. Schröder, and S. Weißmann. Schrödinger’s smoke. ACM Trans. Graph., 35, 7 2016

## Case studies

In this repository, the full code of the case studies investigated by the authors in the supporting papers are available.

| Name                                                                                                    | Reference Paper                                                                                                                                                              |
| ------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Backward Facing Step](https://ermete-lab.github.io/MP-pyISFenix/tutorials/01_ISF_BFS.html)             | [AIP-2024](https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)                                                                  |
| [Flow Over Cylinder](https://ermete-lab.github.io/MP-pyISFenix/tutorials/02_ISF_cyl_2D.html)            | [AIP-2024](https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)                                                                  |
| [Flow Over NACA0012](https://ermete-lab.github.io/MP-pyISFenix/tutorials/03_ISF_NACA0012.html)          | [AIP-2024](https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)                                                                  |
| [Gravity - Tube Bundle](https://ermete-lab.github.io/MP-pyISFenix/tutorials/04_ISF_tube_bundle.html)    | [AIP-2024](https://pubs.aip.org/aip/pof/article/36/1/017138/3132670/A-finite-element-implementation-of-the)                                                                  |
| [Buoyant - Tube Bundle](https://ermete-lab.github.io/MP-pyISFenix/tutorials/05_buoISF_tube_bundle.html) | [UIT-2024](https://www.researchgate.net/publication/381707814_Inclusion_of_the_buoyancy_forces_in_the_Incompressible_Schrodinger_Flow_algorithm_to_simulate_inviscid_fluids) |
| [Nuclear Rod 3D](https://ermete-lab.github.io/MP-pyISFenix/tutorials/06_ISF_NuclearRod3D.html)          | [NUTHOS-2024](https://www.researchgate.net/publication/384146880_Advection-diffusion_of_scalars_with_the_Incompressible_Schrodinger_Flow)                                    |

## How to cite
If you use this code in your research, please cite the following paper

- S. Riva, C. Introini, and A. Cammi, “*A finite element implementation of the incompressible Schroedinger flow method*,” Physics of Fluids, vol. 36, p. 017138, 01 2024


Here in bibtex format

```bibtex
@article{AIP_202401,
    author = {Riva, Stefano and Introini, Carolina and Cammi, Antonio},
    title = "{A finite element implementation of the incompressible Schrödinger flow method}",
    journal = {Physics of Fluids},
    volume = {36},
    number = {1},
    pages = {017138},
    year = {2024},
    month = {01},
    abstract = "{As first proposed by Madelung in 1926, the analogy between quantum mechanics and hydrodynamics has been known for a long time; however, its potentialities and the possibility of using the characteristic equations of quantum mechanics to simulate the behavior of inviscid fluids have not been thoroughly investigated in the past. In this methodology, the incompressible Euler equations are thus substituted by the Schrödinger equation, turning a quasi-linear Partial Differential Equation into a linear one, an algorithm known in the literature as Incompressible Schrödinger Flow. Previous works on the subject used the Fast Fourier Transform method to solve this problem, obtaining promising results, especially in predicting vortex dynamics; this paper aims to implement this novel approach into a Finite Element framework to find a more general formulation better suited for future application on complex geometries and on test cases closer to real-world applications. Simple case studies are presented in this work to analyze the potentialities of this method: the results obtained confirm that this method could potentially have some advantages over traditional Computational Fluid Dynamics method, especially for what concerns computational savings related to the required time discretization, whilst also introducing new aspects of the algorithm, mainly related to boundary conditions, not addressed in previous works.}",
    issn = {1070-6631},
    doi = {10.1063/5.0180356},
    url = {https://doi.org/10.1063/5.0180356},
    eprint = {https://pubs.aip.org/aip/pof/article-pdf/doi/10.1063/5.0180356/18930833/017138\_1\_5.0180356.pdf},
}
```

**Recent works using pyISFenix have extended the code to support advection-diffusion of scalars and buoyancy-like forces at the quantum level**:

- S. Riva, C. Introini, L. Marocco, L. Savoldi, and A. Cammi, “Inclusion of the buoyancy forces in the Incompressible Schroedinger Flow algorithm to simulate inviscid fluids,” in 41st UIT International Heat Transfer Conference, (Naples, Italy), 21 June 2024.
- S. Riva, C. Introini, X. Wang, and A. Cammi, “Advection-Diffusion of Scalars with the Incompressible Schroedinger Flow,” in The 14th International Topical Meeting on Nuclear Reactor Thermal-Hydraulics, Operation, and Safety (NUTHOS24), (Vancouver, USA), August 2024.


## Contact Information

If interested, please contact stefano.riva@polimi.it
