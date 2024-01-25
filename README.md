# MP-pyISFenix

**PYthon framework for the Incompressible Schrodinger Flow using FENIcsX**

**Authors:** Stefano Riva, Carolina Introini, Antonio Cammi

This repository collects some codes implemented to simulate inviscid fluids by means of the ISF technique, firstly developed by [Albert Chern](https://cseweb.ucsd.edu/~alchern/), adopting Finite Element Methods exploiting the [FEniCSx Library (v. 0.6.0)](https://fenicsproject.org/) for Python, extension to v 0.7.2 is under development.

A **webpage** has been created at the following link [https://ermete-lab.github.io/MP-pyISFenix/intro.html](https://ermete-lab.github.io/MP-pyISFenix/intro.html), including the API documentation of the solvers and some tutorials resembling the case studies of the published papers.

--------------------------------

The ISF algorithm is based on the Madelung transform, an analogy between the wave functions governed by the Schrodinger equation and the velocity field for inviscid fluids.

Here's some useful background references:
- E. Madelung, Eine anschauliche deutung der gleichung von schr ̈odinger, Naturwissenschaften 14 (1926) 1004. doi:10.1007/BF01504657.
- A. Chern. Fluid Dynamics with Incompressible Schroedinger Flow. PhD thesis, 2017
- A. Chern, F. Knöppel, U. Pinkall, P. Schröder, and S. Weißmann. Schrödinger’s smoke. ACM Trans. Graph., 35, 7 2016


## How to cite
If you use this code in your research, please cite the following paper

- S. Riva, C. Introini, and A. Cammi, “*A finite element implementation of the incompressible Schroedinger flow method*,” Physics of Fluids, vol. 36, p. 017138, 01 2024


Here in bibtex format

```bibtex
@article{AIP_202401,
    author = {\textbf{\underline{Stefano Riva}} and Introini, Carolina and Cammi, Antonio},
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

## Contact Information

If interested, please contact stefano.riva@polimi.it