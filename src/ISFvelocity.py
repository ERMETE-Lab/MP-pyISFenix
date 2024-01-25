# Schrodinger Velocity
# Author: Stefano Riva, PhD Student, NRG, Politecnico di Milano
# Latest Code Update: 19 January 2023
# Latest Doc  Update: 19 January 2023

import numpy as np
import dolfinx
from dolfinx import fem
from dolfinx.fem import Function, FunctionSpace, dirichletbc, locate_dofs_topological, form, assemble_scalar
import ufl
from ufl import grad, div, nabla_grad, inner, dot, real, FacetNormal, conj
from petsc4py import PETSc

class computeVelocity():
    def __init__(self, domain: dolfinx.mesh.Mesh, degree_psi = 1):
        """
    
        This class implements the calculation of the velocity from a couple of wave functions :math:`[\psi_1, \psi_2]`.
        
        Parameters
        ----------
        domain : dolfinx.mesh.Mesh
            Mesh imported.
        degree_psi : int (Default = 1)
           Degree of the polynomials used for the definition of the functional space.
    
        """
        # Domain
        self.domain = domain
        self.gdim   = self.domain.geometry.dim
        self.fdim   = self.gdim - 1  

        # Definition of the integral differential (both in and on the domain)
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.domain, metadata=metadata)

        # Functional Spaces for velocity and wave function
        self.P1vec = ufl.VectorElement("Lagrange", self.domain.ufl_cell(), 1)
        self.W = FunctionSpace(self.domain, self.P1vec)

        self.P2 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), degree_psi)
        self.V = FunctionSpace(self.domain, self.P2)

        # Definition of the trial and test functions for the velocity
        self.u = ufl.TrialFunction(self.W)
        self.v = ufl.TestFunction(self.W)
        
        # Initialisation of the wave functions
        self.psi_1 = Function(self.V).copy()
        self.psi_2 = Function(self.V).copy()

    def assemble(self, hbar: float, direct : bool = False):
        """
        
        This function assembles the forms of the weak formulations and creates the linear structures for the solution (matrix, rhs and solver).

        The velocity is computed as follows:
        
        .. math::
                \mathbf{u} = \hbar\,\Re\{-(\Psi^*)^T\,i\,\\nabla\Psi\} 
                           = \hbar\,\Re\{-i(\psi_1^*\,\\nabla\psi_1+\psi_2^*\,\\nabla\psi_2)\} 

        Parameters
        ----------
        hbar : float
            Value of :math:`\hbar`.
        direct : boolean (Default: False)
            Boolean variable to identify if a direct solver must be used.
        """
        
        # Assigning hbar
        self.hbar = hbar

        # Defining lhs form
        self.a = form( inner(self.u, self.v) * self.dx)
        
        # Defining rhs for,
        self.L = form( inner( self.hbar * real( - 1.0 * 1j * (conj(self.psi_1) *  grad(self.psi_1) + conj(self.psi_2) *  grad(self.psi_2))),
                              self.v) * self.dx)
        
        # Creating and assembling matrix
        self.A = fem.petsc.assemble_matrix(self.a)
        self.A.assemble()
        
        # Create vector
        self.b = fem.petsc.create_vector(self.L)

        # Define solver 
        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)
        else:
            self.solver.setType(PETSc.KSP.Type.CG)
            self.solver.getPC().setType(PETSc.PC.Type.SOR)

        # Define solution function
        self.solution = Function(self.W).copy()
        
    def solve(self, psi_1: dolfinx.fem.Function, psi_2: dolfinx.fem.Function):
        """
        
        This function computes the velocity given a couple of wave functions.
        
        Parameters
        ----------
        psi_1 : dolfinx.fem.Function
            First component of the wave function :math:`\psi_1`.
        psi_2 : dolfinx.fem.Function
            Second component of the wave function :math:`\psi_2`.
            
        """
        
        # Assign the wave function to the structures of the class
        with psi_1.vector.localForm() as loc_, self.psi_1.vector.localForm() as loc_n:
            loc_.copy(loc_n)
        with psi_2.vector.localForm() as loc_, self.psi_2.vector.localForm() as loc_n:
            loc_.copy(loc_n)

        # Assembling the rhs vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)

        # Solving the linear system
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()

        return self.solution