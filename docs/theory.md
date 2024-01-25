# Theory

In this page, some theoretical aspects are discussed regarding the idea behind Madelung transformation, the ISF algorithm and some elements on functional analysis and weak formulation.
For a complete description of the mathematics of the methods please refer to the ISF papers {cite:p}`Riva_ISF, Chern2016, Chern2017` or to numerical analysis books {cite:p}`Gazzola_A3, quarteroni_2016, QuarteroniMateNumerica`. For the Finite Element implementation and the derivation of the weak formulations for this problem see the paper {cite:p}`Riva_ISF`.

## Madelung Transformation

The *Incompressible Schrodinger Flow* is a novel numerical technique, developed in {cite}`Chern2016` and {cite}`Chern2017`, able to describe inviscid fluids behaviour using the Schrodinger equation. This section gives some basic concepts on the analogy between hydrodynamics and quantum mechanics {cite}`Madelung1926`, onto which ISF is based, along with the description of the solution algorithm.

The Schrodinger equation is able to describe the time evolution of this physical state, subjected to a certain potential field $p$
\begin{equation*}
    i\hbar \frac{\partial\psi}{\partial t} = -\frac{\hbar^2}{2}\Delta \psi + p \psi.
\end{equation*}
In the context of ISF, the wave function is used to describe inviscid flows, the probabilistic interpretation of the quantum mechanics is no more used, instead the ISF is set in the analogy between hydrodynamics and quantum mechanics is no more used, instead the ISF is set in the analogy between hydrodynamics and quantum mechanics, proposed by Madelung in 1926. The state of the system can be described by a 2 components wave function {cite:p}`Chern2016`, i.e. $\Psi=[\psi_1,\psi_2]^T\in\mathbb{C}^2$, which is the solution of the following Schrodinger equation, i.e.
\begin{equation*}
    i\hbar\frac{\partial\Psi}{\partial t} =- \frac{\hbar^2}{2}\Delta \Psi+p\Psi\longleftrightarrow
    i\hbar\frac{\partial}{\partial t}\left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right] = -\frac{\hbar^2}{2}\Delta \left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right]+p\left[\begin{array}{cc} \psi_1 \\ \psi_2\end{array}\right].
\end{equation*}

```{note}
If the flow field were described by a single wave function wave function $\psi\in\mathbb{C}$, the associated velocity field would be characterised by a null vorticity, which is generally not true.
```

The wave function can be translated into a flow field as
\begin{equation*}
    \mathbf{u} = \hbar\mbox{Re}\left\{-\left(\Psi^T\right)^*i\,\nabla\Psi\right\}.
\end{equation*}

It can be shown that the real and the imaginary part of this equation are linked to the continuity and momentum equation of the incompressible Euler equations
\begin{equation*}
    \left\{
    \begin{array}{l}
        \nabla\cdot \mathbf{u} = 0 \\
        \displaystyle\frac{\partial\mathbf{u}}{\partial t}+\left(\mathbf{u}\cdot \nabla\right)\mathbf{u} = -\nabla p
    \end{array}
    \right.
\end{equation*}

The viscosity, thus the dissipation, is not considered and the only existing external force is given by the potential $-\nabla p$. The spatial operators on the wave function $\Psi$, i.e. $-\frac{\hbar^2}{2}\Delta$ and $p$, have a physical interpretation in the Euler equations on $\mathbf{u}$, i.e. the advection $\left(\mathbf{u}\cdot \nabla\right)\mathbf{u}$ and the potential forces $-\nabla p$.

In the end, the incompressibility constraint $\nabla \cdot \mathbf{u}=0$ can be equivalently imposed to the wave function as
\begin{equation*}
    \mbox{Re}\left\{\left(\Psi^T\right)^*i\,\Delta\Psi\right\} = 0.
    \label{eqn:SchrodingerIncCoinstraint}
\end{equation*}

## ISF algorithm

Let $\Omega$ be the spatial domain, $\partial \Omega$ its boundary, usually composed by three different components $\Gamma_{in}\cup \Gamma_w\cup \Gamma_o$ and let $\mathcal{T} = [0, T]$ be the time interval considered. Given an initial condition $\Psi^0$, the Schrodinger equation with the incompressibility constraint, can be solved with a first order time splitting method

```{note}
A similar solution strategy is also employed in CFD, a very important example is the pioneering work of {cite}`Chorin1968`.
```

At each time step $t_n\rightarrow t_{n+1}$, the algorithm, in the space-continuous case solve the following problems:

- Schrodinger problem
- Normalisation
- Poisson equation (*pressure projection*)
- Phase Shift and velocity update

### Schrodinger problem
The prediction step consists of a free-particle Schrodinger equation
\begin{equation*}
    \frac{\tilde{\Psi}-\Psi^n}{\Delta t} = \frac{i \hbar}{2}\Delta \tilde{\Psi}
\end{equation*}
with the associated boundary conditions (a plane wave at the inlet and homogeneous Neumann boundaries for the others), i.e. 
\begin{equation*}
     \left\{
        \begin{array}{ll}
            \displaystyle\tilde{\Psi} = e^{i\left(\mathbf{k}\cdot \mathbf{x}-\omega t\right)}\cdot \left[c_1, c_2\right]^T & \text{ on }\Gamma_{in}\\
            \nabla\tilde{\Psi}\cdot \mathbf{n}=0 & \text{ on }\Gamma_w\cup \Gamma_o
        \end{array}
        \right.
    \end{equation*}
with $c_1,c_2 \in\mathbb{C}$ such that $\sqrt{|c_1|^2+|c_2|^2}=1$.

### Normalisation
The wave function needs also to be normalised to 1, since this property is not necessarily conserved in the time-discrete case, i.e.
\begin{equation*}
    \tilde{\Psi} \longleftarrow\frac{\tilde{\Psi}}{\;||\tilde{\Psi}||_2} =
    \frac{\tilde{\Psi}}{\sqrt{\tilde{\psi}_1^*\tilde{\psi}_1+\tilde{\psi}_2^*\tilde{\psi}_2}}  .
\end{equation*}

### Pressure Projection
The correction phase involves a Poisson problem to enforce the incompressibility constraint. The unknown is $\varphi$, which is referred to as {pressure}, even though its units of measure are $[m^2/s]$. This quantity has been introduced following the idea developed in {cite}`Chern2017`, which provides a complete discussion on this choice. Thus, the following Poisson problem is solved
\begin{equation*}
        \left\{
    \begin{array}{ll}
       \Delta \varphi = \nabla \cdot \tilde{\mathbf{u}} & \mbox{ in }\Omega \\
       \varphi = 0 & \mbox{ on } \Gamma_{in}\\
        \nabla\varphi\cdot \mathbf{n}=0 & \mbox{ on } \Gamma_w\\ 
        \text{Outlet BC} & \mbox{ on } \Gamma_o
    \end{array}
    \right. 
    \label{eqn: Poisson-Strong}
\end{equation*}
in which $\tilde{\mathbf{u}}$ is the velocity field computed from $\tilde{\Psi}$. The outlet boundary condition can be either a non-homogeneous Neumann or an homogeneous Dirichlet, refer to {cite:p}`Riva_ISF` for the details.

```{note}
This idea is very similar to the standard CFD splitting method, in which a Poisson problem is solved to enforce the incompressibility.
```

### Phase Shift and Velocity Update

In the end, the wave function and the velocity are updated, i.e.
\begin{equation*}
\begin{split}
    \Psi^{n+1} &= e^{-i\frac{\varphi}{\hbar}}\tilde{\Psi}\\
    \mathbf{u}^{n+1}&= \hbar\,\Re\{-(\Psi^{n+1,*})^T\,i\,\nabla\Psi^{n+1}\} \\
                    &= \hbar\,\Re\{-i(\psi_1^{n+1,*}\,\nabla\psi_1^{n+1}+\psi_2^{n+1,*}\,\nabla\psi_2^{n+1})\} 
\end{split}
\end{equation*}

## Elements of Functional Analysis

Before entering into the details of the weak formulations, some definitions must be presented. First, let $L^2(\Omega)$ be an Hilbert space such that
\begin{equation*}
   u\in L^2(\Omega) \Longleftrightarrow \int_\Omega |u|^2\,d\Omega <\infty;
\end{equation*}
endowed with a scalar product
\begin{equation*}
    \langle{u,v}\rangle = \int_\Omega v^*u\,d\Omega \qquad\qquad \forall u,v\in L^2(\Omega).
\end{equation*}
The *weak* solutions are usually searched in the Sobolev functional spaces $\mathcal{H}^p(\Omega)$ {cite:p}`salsa, Gazzola_A3` such that
\begin{equation*}
    \mathcal{H}^p(\Omega) = \left\{u\in L^2(\Omega)\,:\,\int_\Omega |{D^{\boldsymbol{\alpha}}u}|^2\,d\Omega <\infty \right\},
\end{equation*}
given $D^{\boldsymbol{\alpha}}$ the multi-index derivative {cite:p}`Gazzola_A3` of order $p$. Within this mathematical framework, the *weak* formulation can be obtained by testing the *strong* form against a {test function} living in a suitable Sobolev space. In the following, the weak formulations for the Schrodinger and Poisson problems will be derived.

### Useful formulas

Let $\mathbf{u}\in[\mathcal{H}^1(\Omega)]^d$ a vector function, the Gauss (divergence) theorem states:

$$
\int_\Omega \nabla \cdot \mathbf{u}\,d\Omega = \int_{\partial\Omega}\mathbf{u}\cdot \mathbf{n}\,d\sigma
$$

From this theorem, the following formulas can be derived
\begin{equation*}
\int_\Omega (\nabla \cdot \mathbf{u}) \,p\,d\Omega = 
- \int_\Omega \mathbf{u} \cdot \nabla p\,d\Omega + \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
\end{equation*}
given $p\in\mathcal{H}^1(\Omega)$.

````{prf:proof}
Recalling that $ \nabla\cdot (p\,\mathbf{u})=(\nabla \cdot \mathbf{u}) \,p + \mathbf{u} \cdot \nabla p$ we can rewrite the formula as

$$
\int_\Omega \left[(\nabla \cdot \mathbf{u}) \,p + \mathbf{u} \cdot \nabla p\right] \,d\Omega = 
\int_\Omega \nabla \cdot (p\mathbf{u})\,d\Omega = \int_{\partial\Omega}(\mathbf{u}\cdot \mathbf{n})\,p\,d\sigma
$$
from which thesis follows.
````

## Weak Formulations

Let $\Omega\subset\mathbb{R}^{d}$ with $d=2,3$ and $\partial \Omega = \Gamma_D\cup\Gamma_N$ its boundary, let us consider the following problem (in strong form):
\begin{equation*}
\left\{
    \begin{array}{ll}
        -\Delta u = f & \mbox{in } \Omega\\
        u = u_D & \mbox{on } \Gamma_D\\
        \displaystyle \frac{\partial u}{\partial \mathbf{n}} = g_N & \mbox{on } \Gamma_N
    \end{array}
\right.
\end{equation*}
We have a PDE to solve with its boundary condition (BC) associated. The strong solution must be differentiable twice in classical sense and, in general, it is not always possible to prove the well-posedness of these problems. Therefore, it comes the necessity to look for the solution into a broader functional space. Without entering into too much details on this topic, this sections aims to present the typical procedure to be followed when a weak formulation must be derived (its necessary to assemble the algebraic system coming from the Finite element method).

Let $\mathcal{V}\subset\mathcal{H}^1, \mathcal{V}_0\subset\mathcal{H}^1$ the trial and test space defined as

$$
\mathcal{V} = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = u_D\right\}\qquad 
\mathcal{V}_0 = \left\{v\in\mathcal{H}^1:\;\left. v\right|_{\Gamma_D} = 0\right\}
$$

Given $v\in\mathcal{V}_0$, the PDE is multipied (in $L^2$ sense) by $v$

$$
\int_\Omega -\Delta u\,v \,d\Omega=\int_\Omega f\,v\,d\Omega\qquad v\in\mathcal{V}_0
$$
and by applying the integration by parts formula and the effect of the BCs
\begin{equation*}
\int_\Omega \nabla u\cdot \nabla v \,d\Omega=\int_\Omega f\,v\,d\Omega+\int_{\Gamma_N}g\,v\,d\sigma\qquad v\in\mathcal{V}_0
\end{equation*}

### Galerkin approximation
Let us introduce a computational grid $\mathcal{T}_h$ (dim$\mathcal{T}_h = \mathcal{N}_h$) of the domain $\Omega$ and the finite dimensional spaces $\mathcal{V}_h\subset\mathcal{V}$ and $\mathcal{V}_{h,0}\subset\mathcal{V}_0$. The Galerkin approximation {cite}`quarteroni_2016, QuarteroniMateNumerica` reads: find $u_h\in\mathcal{V}_h$ s.t.
\begin{equation*}
\int_\Omega \nabla u_h\cdot \nabla v_h \,d\Omega=\int_\Omega f_h\,v_h\,d\Omega+\int_{\Gamma_N}g\,v_h\,d\sigma\qquad v_h\in\mathcal{V}_{h,0}
\end{equation*}

### Finite element algebraic formulation
The finite dimensional formulation is necessary to assemble the associated linear system. 
Let us consider the homogeneous Dirichlet case, $u_D = 0$ (to which we can
always reduce through a lifting procedure), and let us introduce a basis $\left\{\phi_j\right\}_{j=1}^{\mathcal{N}_h}$ such that any function $u_h\in\mathcal{V}_h$ can be expressed as a linear combination of the basis functions

$$
u_h(\mathbf{x}) = \sum_{j=1}^{\mathcal{N}_h}\alpha_j \phi_j(\mathbf{x})
$$

```{note}

The choice of the basis functions results in different versions of the algebraic system. The most natural choice is given by the **hat functions** {cite:p}`QuarteroniMateNumerica`.
```

By taking $v_h=\phi_k$ as the test function, the Galerkin problem reduces to
\begin{equation*}
\int_\Omega \nabla u_h\cdot \nabla \phi_k \,d\Omega=\int_\Omega f_h\,\phi_k\,d\Omega+\int_{\Gamma_N}g\,\phi_k\,d\sigma\qquad k = 1, \dots, \mathcal{N}_h
\end{equation*}
then, the linear expansion can be introduced into the problem which results

$$
\sum_{j=1}^{\mathcal{N}_h}\alpha_j \int_\Omega \nabla \phi_j\cdot \nabla \phi_k \,d\Omega=\int_\Omega f_h\,\phi_k\,d\Omega+\int_{\Gamma_N}g\,\phi_k\,d\sigma\qquad k = 1, \dots, \mathcal{N}_h
$$
By defining the matrix $A_{kj} = \displaystyle\int_\Omega \nabla \phi_j\cdot \nabla \phi_k \,d\Omega$ and the vectors $f_k =  \displaystyle\int_\Omega f_h\,\phi_k\,d\Omega$ and $g_k =  \displaystyle\int_{\Gamma_N}g\,\phi_k\,d\sigma$, the problem into the following linear system in the unknowns $\alpha_j$
\begin{equation*}
\sum_{j=1}^{\mathcal{N}_h}A_{kj}\alpha_j = f_k + g_k\qquad k = 1, \dots, \mathcal{N}_h
\end{equation*}