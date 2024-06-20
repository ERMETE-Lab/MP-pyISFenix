# Schrodinger Advection Step
# Author: Stefano Riva, PhD Student, NRG, Politecnico di Milano
# Latest Code Update: 24 November 2023
# Latest Doc  Update: 24 November 2023

from math import isclose
import numpy as np
import dolfinx
from dolfinx import fem
from dolfinx.fem import petsc, Function, FunctionSpace, dirichletbc, locate_dofs_topological, form, assemble_scalar
import ufl
from ufl import grad, div, nabla_grad, inner, dot, real, FacetNormal
from petsc4py import PETSc

class PlaneWave():
    """
    A class to generate a plane wave for the wave function, defined as
    
    .. math::
            \psi(\mathbf{x}, t) = c \cdot e^{ i\,(\mathbf{k} \cdot \mathbf{x} - \omega \, t) }
    
    in which :math:`c\in\mathbb{C}` is a normalisation constant, :math:`\mathbf{k}` is the wave vector and :math:`\omega` is the frequency, defined as
    
    .. math::
        \mathbf{k} = h^{-1}\cdot \mathbf{u} \qquad \qquad \omega = (\mathbf{u}\cdot \mathbf{u}) / (2\hbar)
        
        
    The plane wave is used to assign a fixed velocity field to a certain region.
    
    Parameters
    ----------
    t : float
      Initial Time.
    kVec : np.ndarray
      Wave vector.
    omega : float
      Frequency of the wave.
    hbar : float
      Value of :math:`\hbar`.
    c : complex
      Normalisation constant.

    """
    def __init__(self, t: float, kVec: np.ndarray, omega: float, hbar: float, c: complex):
        self.t = t
        self.hbar = hbar
        self.kVec = kVec
        self.omega = omega
        self.c = c
    def __call__(self, x):
        """
        Defines the class as callable, this is required to interpolate this into a `dolfinx.fem.Function` object.
        """
        
        # Definition of the spatial contribution of the plane wave
        spatial = sum([self.kVec[ii] * x[ii] for ii in range(len(self.kVec))])
        
        return self.c * np.exp( 1j * ( spatial - self.omega * self.t))

class schrodinger():
    """
    
    This class implements the prediction step of the Incompressible Schrodinger Flow, in particular the free-particle Schrodinger equation is solved.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
      Mesh imported.
    ft : dolfinx.cpp.mesh.MeshTags_int32
      Face tags of the boundaries.
    boundary_marks : dict
      Dictionary with the markers for the boundaries.
    degree_psi : int (Default = 1)
      Degree of the polynomials used for the definition of the functional space.
    c : np.ndarray (Default = :math:`[1, 0.01]`)
      Constant, in the `__init__` is normalised such that :math:`\sqrt{c_1^*c_1+c_2^*c_2}`.
    """
    def __init__(self, domain: dolfinx.mesh.Mesh, ft: dolfinx.cpp.mesh.MeshTags_int32, boundary_marks: dict, 
                 degree_psi: int = 1, c: np.ndarray = np.array([1.0 + 0 * 1j, 0.01 + 0 * 1j])):
        
        # Domain
        self.domain = domain
        self.ft = ft
        self.gdim   = self.domain.geometry.dim
        self.fdim   = self.gdim - 1  

        # Definition of the integral differential (both in and on the domain)
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.domain, metadata=metadata)
        self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=self.ft, metadata=metadata)

        # Assign markers: only inlet is needed, on the other boundaries homogeneous Neumann is imposed
        self.inl_marker  = boundary_marks['inlet']

        # Functional Spaces for the wave function
        self.P2 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), degree_psi)
        self.mixEl = ufl.MixedElement(self.P2, self.P2)
        self.V = FunctionSpace(self.domain, self.mixEl)
        
        # Functional Spaces for the temperature - buoyancy effects
        self.Qn = FunctionSpace(self.domain, ("Lagrange", 1))

        # Definition of the trial and test functions for the wave function
        (self.psi1, self.psi2) = ufl.TrialFunctions(self.V)
        (self.v1,   self.v2)   = ufl.TestFunctions(self.V)

        # Initialisation
        self.t = 0.
        
        self.c1 = c[0]
        self.c2 = c[1]

        constantNorm = np.sqrt(np.conj(self.c1) * self.c1 + np.conj(self.c2) * self.c2)
        self.c1 /= constantNorm
        self.c2 /= constantNorm
        
        self.psi1_old = Function(self.V.sub(0).collapse()[0]).copy()
        self.psi2_old = Function(self.V.sub(1).collapse()[0]).copy()
        
        self.psi1_old.interpolate(lambda x: self.c1 + x[0] * 0.0 + x[1] * 0.0)
        self.psi2_old.interpolate(lambda x: self.c2 + x[0] * 0.0 + x[1] * 0.0) 
        # self.psi1_old.x.set(self.c1)
        # self.psi2_old.x.set(self.c2) 

    def updateInlet(self):
        """
        This function updates the inlet condition of the plane wave.
        """
        
        self.psi1_in.interpolate(self.wave1)
        self.psi2_in.interpolate(self.wave2)

    def assemble(self,  phys_parameters: dict, dt: float = 1e-2, direct : bool = False,
                        gravity_dict: dict = {'gravity': np.array([0, 0, 0]), 'beta': 0, 'Tref': 0} ):
        """
        
        This function assign the boundary condition (if inlet is present), assembles the forms of the weak formulations and creates the linear structures for the solution (matrix, rhs and solver).

        Parameters
        ----------
        phys_parameters : dict
            Dictionary with the physical parameters.
        dt : float (Default: 1e-2)
            Time step.
        direct : boolean (Default: False)
            Boolean variable to identify if a direct solver must be used.
        gravity_dict : dict (Default: None)
            Dictionary with the gravity parameters.
        """
        
        # Assign the time step
        self.dt = dt
        
        # Assign the hbar value
        self.hbar = phys_parameters['hbar']
        
        # Assign the inlet BC, if required
        if self.inl_marker is not None:
            
            # Defining wave vector and frequency
            self.kVec = phys_parameters['uin'] / self.hbar
            self.omega = np.dot(phys_parameters['uin'], phys_parameters['uin']) / ( 2 * self.hbar)

            # Defining inlet waves
            self.psi1_in = Function(self.V.sub(0).collapse()[0]).copy()
            self.psi2_in = Function(self.V.sub(1).collapse()[0]).copy()
            self.wave1 = PlaneWave(self.t, self.kVec, self.omega, self.hbar, self.c1)
            self.wave2 = PlaneWave(self.t, self.kVec, self.omega, self.hbar, self.c2)
            
            # Updating value of the inlet waves
            self.updateInlet()
            
            # Define the inlet boundary conditions
            self.bc_in1 = dirichletbc(self.psi1_in, locate_dofs_topological((self.V.sub(0), self.V.sub(0).collapse()[0]), self.fdim, self.ft.find(self.inl_marker)), self.V.sub(0))
            self.bc_in2 = dirichletbc(self.psi2_in, locate_dofs_topological((self.V.sub(1), self.V.sub(1).collapse()[0]), self.fdim, self.ft.find(self.inl_marker)), self.V.sub(1))

            self.bcs = [self.bc_in1, self.bc_in2]
        else:
            # If there is no inlet the whole boundary has homogeneous Neumann imposed.
            self.bcs = []

        # Define the left-hand side of the form
        self.lhs  = 1. / self.dt * inner(self.psi1, self.v1) * self.dx
        self.lhs += 1. / self.dt * inner(self.psi2, self.v2) * self.dx
        self.lhs += 1j * self.hbar / 2. * inner(grad(self.psi1), grad(self.v1)) * self.dx
        self.lhs += 1j * self.hbar / 2. * inner(grad(self.psi2), grad(self.v2)) * self.dx

        ## If the gravitational field is considered an additional potential term is added
        if sum(abs(gravity_dict['gravity'])) > 0:
            
            # Defining gravity functions
            self.g1 = Function(self.V.sub(0).collapse()[0]).copy()
            self.g2 = Function(self.V.sub(1).collapse()[0]).copy()
            
            self.g1.interpolate(lambda xx: +(xx[0] * gravity_dict['gravity'][0] + xx[1] * gravity_dict['gravity'][1] + xx[2] * gravity_dict['gravity'][2]))
            self.g2.interpolate(lambda xx: -(xx[0] * gravity_dict['gravity'][0] + xx[1] * gravity_dict['gravity'][1] + xx[2] * gravity_dict['gravity'][2]))
            
            # Defining the temperature functions
            self.Tref = Function(self.Qn).copy()
            self.Told = Function(self.Qn).copy()
            self.Tref.interpolate(lambda x: gravity_dict['Tref'] + x[0] * 0.0 + x[1] * 0.0)
            self.Told.interpolate(lambda x: gravity_dict['Tref'] + x[0] * 0.0 + x[1] * 0.0)
            
            self.beta = PETSc.ScalarType(gravity_dict['beta'])
        
            # Defining the buoyancy term and add to the lhs
            self._g1 = self.g1 * ( 1 - self.beta * (self.Told - self.Tref) )
            self._g2 = self.g2 * ( 1 - self.beta * (self.Told - self.Tref) )
        
            self.lhs += inner(self._g1 * self.psi1, self.v1) * self.dx
            self.lhs += inner(self._g2 * self.psi2, self.v2) * self.dx

        # Forming the lhs
        self.a = form( self.lhs )

        # Defining the right-hand side of the form 
        self.L = form( 1. / self.dt * inner(self.psi1_old, self.v1) * self.dx +
                       1. / self.dt * inner(self.psi2_old, self.v2) * self.dx)

        # Creating matrix and vector
        self.A = petsc.create_matrix(self.a)
        petsc.assemble_matrix(self.A, self.a, self.bcs)
        self.A.assemble()
        self.b = fem.petsc.create_vector(self.L)

        # Define solver 
        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)
        else:
            self.solver.setType(PETSc.KSP.Type.BCGS)
            self.solver.getPC().setType(PETSc.PC.Type.JACOBI)
        
        # Define solution structures
        self.solution = Function(self.V).copy()
        # self.psi1_norm, self.psi2_norm = self.solution.sub(0).collapse(), self.solution.sub(1).collapse()
        
    def advect(self, temperature: dolfinx.fem.Function = None):
        """
        
        This function advances in time by solving the linear problem at each time step.
        
        Parameters
        ----------
        temperature : dolfinx.fem.Function (Default = None)
            Temperature at time :math:`t_n`, if None there is no update.
            
        """
        
        self.t += self.dt
        
        # Updating BC
        if self.inl_marker is not None:
          self.wave1.t = self.t
          self.wave2.t = self.t
          self.updateInlet()

        # Updating Temperature - If not None
        if temperature is not None:
            self.Told.x.array[:] = temperature.x.array
            self.A.zeroEntries()
            fem.petsc.assemble_matrix(self.A, self.a, self.bcs)
            self.A.assemble()
        
        # Assembling right-hand side vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        fem.petsc.apply_lifting(self.b, [self.a], [self.bcs])
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(self.b, self.bcs)

        # Solving the linear problem and returning solution
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()

        return self.solution.sub(0).collapse(), self.solution.sub(1).collapse()
      
##########################################################################################################################################
###################                                Additional routines                                                    ################
##########################################################################################################################################

#     def extractInlet(self, wave_fun, in_size):
    
#         Ny = 50
#         y_grid = np.linspace(0., in_size, Ny)

#         points = np.zeros((3, Ny))
#         points[0, :] = 0.
#         points[1, :] = y_grid

#         bb_tree = dolfinx.geometry.BoundingBoxTree(self.domain, self.domain.topology.dim)
#         cells = []
#         points_on_proc = []
#         cell_candidates = dolfinx.geometry.compute_collisions(bb_tree, points.T)
#         colliding_cells = dolfinx.geometry.compute_colliding_cells(self.domain, cell_candidates, points.T)
#         for i, point in enumerate(points.T):
#             if len(colliding_cells.links(i))>0:
#                 points_on_proc.append(point)
#                 cells.append(colliding_cells.links(i)[0])
#         xPlot = np.array(points_on_proc, dtype=np.float64)

#         in_wave = wave_fun.eval(xPlot, cells).flatten()
#         return in_wave[int(Ny/2)]
    
# class normalisation():
#     """_summary_
#     """
#     def __init__(self, domain, ft, boundary_marks : dict, degree_psi = 2):
        
#         # Domain
#         self.domain = domain
#         self.ft = ft
#         self.gdim   = self.domain.geometry.dim
#         self.fdim   = self.gdim - 1  

#         metadata = {"quadrature_degree": 4}
#         self.dx = ufl.Measure("dx", domain=self.domain, metadata=metadata)
#         self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=self.ft, metadata=metadata)

#         # Assign markers
#         self.inl_marker  = boundary_marks['inlet']
#         self.out_marker  = boundary_marks['outlet']

#         # Functional Spaces
#         self.P2 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), degree_psi)
#         self.mixEl = ufl.MixedElement(self.P2, self.P2)
#         self.V = FunctionSpace(self.domain, self.mixEl)

#         (self.psi1, self.psi2) = ufl.TrialFunctions(self.V)
#         (self.v1,   self.v2)   = ufl.TestFunctions(self.V)

#     def assemble(self):
        
#         self.psiTilde = Function(self.V)
#         (self.psiTilde_1, self.psiTilde_2) = (self.psiTilde.sub(0).collapse(), self.psiTilde.sub(1).collapse())

#         self.a = form( inner(self.psi1, self.v1) * self.dx + inner(self.psi2, self.v2) * self.dx )
#         self.L = form( inner(self.psiTilde_1 * ufl.algebra.Power( inner(self.psiTilde_1, self.psiTilde_1) + 
#                                                                   inner(self.psiTilde_2, self.psiTilde_2), -0.5),
#                              self.v1) * self.dx +
#                        inner(self.psiTilde_2 * ufl.algebra.Power( inner(self.psiTilde_1, self.psiTilde_1) + 
#                                                                   inner(self.psiTilde_2, self.psiTilde_2), -0.5),
#                              self.v2) * self.dx )

#         self.A = fem.petsc.assemble_matrix(self.a)
#         self.A.assemble()
#         self.b = fem.petsc.create_vector(self.L)

#         self.solver = PETSc.KSP().create(self.domain.comm)
#         self.solver.setOperators(self.A)
#         self.solver.setType(PETSc.KSP.Type.CG)
#         self.solver.getPC().setType(PETSc.PC.Type.SOR)
#         # self.solver.setType(PETSc.KSP.Type.PREONLY)
#         # self.solver.getPC().setType(PETSc.PC.Type.LU)

#     def solve(self, psiTilde_1, psiTilde_2):

#         with psiTilde_1.vector.localForm() as loc_, self.psiTilde_1.vector.localForm() as loc_n:
#             loc_.copy(loc_n)
#         with psiTilde_2.vector.localForm() as loc_, self.psiTilde_2.vector.localForm() as loc_n:
#             loc_.copy(loc_n)
        
#         with self.b.localForm() as loc:
#             loc.set(0)
#         fem.petsc.assemble_vector(self.b, self.L)
#         self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)

#         solution = Function(self.V).copy()
#         self.solver.solve(self.b, solution.vector)

#         return (solution.sub(0).collapse(), solution.sub(1).collapse())