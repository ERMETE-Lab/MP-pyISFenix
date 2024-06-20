import numpy as np
import dolfinx
from dolfinx import fem
from dolfinx.fem import Function, FunctionSpace, dirichletbc, locate_dofs_topological, form, assemble_scalar
import ufl
from ufl import grad, div, nabla_grad, inner, dot, real, FacetNormal
from petsc4py import PETSc
from mpi4py import MPI

class poisson():
    """
    This class implements the projection step of the Incompressible Schrodinger Flow, in particular the solution of the Poisson problem.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
      Mesh imported.
    ft : dolfinx.cpp.mesh.MeshTags_int32
      Face tags of the boundaries.
    boundary_marks : dict
      Dictionary with the markers for the boundaries.
    
    """
    def __init__(self, domain: dolfinx.mesh.Mesh, ft: dolfinx.cpp.mesh.MeshTags_int32, boundary_marks: dict):
        
        # Domain
        self.domain = domain
        self.ft = ft
        self.gdim   = self.domain.geometry.dim
        self.fdim   = self.gdim - 1  
        
        # Definition of the integral differential (both in and on the domain)
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.domain, metadata=metadata)
        self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=self.ft, metadata=metadata)

        # Assign markers: only inlet/outlet is needed, on the other boundaries homogeneous Neumann is imposed
        self.inl_marker  = boundary_marks['inlet']
        self.out_marker  = boundary_marks['outlet']

        # Functional Spaces for the pressure and velocity
        self.P1 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), 1)
        self.Q  = FunctionSpace(self.domain, self.P1)

        self.P1vec = ufl.VectorElement("Lagrange", self.domain.ufl_cell(), 1)
        self.W = FunctionSpace(self.domain, self.P1vec)
        
        # Definition of the trial and test functions for the pressure
        self.p = ufl.TrialFunction(self.Q)
        self.q = ufl.TestFunction(self.Q)
        
        # Initialize velocity field
        self.uTilde = Function(self.W).copy()

    def assemble(self, bc_values: dict, gauss : bool = False, direct : bool = False):
        """
        
        This function assign the boundary condition (if inlet is present), assembles the forms of the weak formulations and creates the linear structures for the solution (matrix, rhs and solver).

        The weak formulation is as follows:
        
        .. math::
            \int_\Omega \\nabla \\varphi \cdot \\nabla q\, d\Omega = - \int_\Omega q\,\\nabla \cdot \\tilde{\mathbf{u}}\, d\Omega

        there is an option to apply the integration by parts on the rhs as
        
        .. math::
            \int_\Omega \\nabla \\varphi \cdot \\nabla q\, d\Omega =  \int_\Omega \\nabla q\cdot \\tilde{\mathbf{u}}\, d\Omega - \int_{\partial \Omega} \\tilde{\mathbf{u}}\cdot \mathbf{n}\,q\,d\sigma
        

        Parameters
        ----------
        bc_values : dict
            Dictionary with the boundary conditions parameters.
        gauss : boolean (Default: False)
            If `True`, the integration by parts is applied also to the right-hand side.
        direct : boolean (Default: False)
            Boolean variable to identify if a direct solver must be used.
        """
        
        # Defining boundary conditions
        if bc_values['fully_Neumann']:
            self.bcs = []
            self.g = None
        else:
            self.inl_p = Function(self.Q).copy()
            # self.inl_p.x.set(0.) 

            self.bc_in  = dirichletbc(self.inl_p, locate_dofs_topological(self.Q, self.fdim, self.ft.find(self.inl_marker)))
            self.bcs = [self.bc_in] 

            if bc_values['is_dirichlet']:
                self.out_p = Function(self.Q).copy()
                # self.out_p.x.set(0.0)
                self.bc_out = dirichletbc(self.out_p, locate_dofs_topological(self.Q, self.fdim, self.ft.find(self.out_marker)))
                self.bcs.append(self.bc_out)
                self.g = None
                
            else:
                
                self.inletMeasure  = assemble_scalar( form(1. * self.ds(self.inl_marker)) )
                self.outletMeasure = assemble_scalar( form(1. * self.ds(self.out_marker)) )
                
                self.g = Function(self.Q).copy()
                self.g.interpolate(lambda x: - self.inletMeasure / self.outletMeasure * bc_values['uin'][0] + 0.0 * x[0])
                # self.g.x.set( - self.inletMeasure / self.outletMeasure * bc_values['uin'][0] ) 

        # Define the lhs of the weak formulation
        self.a = form( inner(grad(self.p), grad(self.q)) * self.dx)
        
        # Define the rhs of the weak formulation - 
        if gauss:
            self.n = FacetNormal(self.domain)
            self.rhs = inner(self.uTilde, grad(self.q)) * self.dx - inner(dot(self.uTilde, self.n), self.q) * self.ds
        else:
            self.rhs = - inner(div(self.uTilde), self.q) * self.dx
        
        # Adding boundary term to the rhs if Neumann is used
        if self.g is not None:
            #  print('Using Neumann')
             self.rhs += inner(self.g, self.q) * self.ds(self.out_marker)
        
        self.L = form( self.rhs )

        
        # Creating and assembling matrix
        self.A = fem.petsc.assemble_matrix(self.a, bcs = self.bcs)
        self.A.assemble()
        
        # Create vector
        self.b = fem.petsc.create_vector(self.L)

        # Set null space if pure Neumann problem
        if bc_values['fully_Neumann']:
            self.nullspace = PETSc.NullSpace().create(constant=True,comm=MPI.COMM_WORLD)
            self.A.setNullSpace(self.nullspace)
            self.nullspace.remove(self.b)
            
        # Define solver 
        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)
        else:
            self.solver.setType(PETSc.KSP.Type.GMRES)
            # self.solver.getPC().setType(PETSc.PC.Type.ILU)
            # self.solver.getPC().setType(PETSc.PC.Type.QR)
            self.solver.getPC().setType(PETSc.PC.Type.SOR)
            # self.solver.getPC().setType(PETSc.PC.Type.JACOBI) # a bit slow 
        
        # Define solution function
        self.solution = Function(self.Q).copy()
        
    def solve(self, uTilde: dolfinx.fem.Function):
        """
        
        This function solves the Poisson problem.
        
        Parameters
        ----------
        uTilde : dolfinx.fem.Function
            Velocity field coming from the prediction step.
            
        """
        
        # Assign the velocity to the structures of the class
        with uTilde.vector.localForm() as loc_, self.uTilde.vector.localForm() as loc_n:
            loc_.copy(loc_n)

        # Assembling the rhs vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        fem.petsc.apply_lifting(self.b, [self.a], [self.bcs])
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(self.b, self.bcs)

        # Solving the linear system
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()

        return self.solution
    
class phaseShift():
    """
    
    This class implements the phase shift of the wave functions :math:`[\psi_1, \psi_2]`.
    
    Parameters
    ----------
    domain : dolfinx.mesh.Mesh
        Mesh imported.
    degree_psi : int (Default = 1)
        Degree of the polynomials used for the definition of the functional space.

    """
    def __init__(self, domain: dolfinx.mesh.Mesh, degree_psi = 1):
        
        # Domain
        self.domain = domain
        self.gdim   = self.domain.geometry.dim
        self.fdim   = self.gdim - 1  

        # Definition of the integral differential (both in and on the domain)
        metadata = {"quadrature_degree": 4}
        self.dx = ufl.Measure("dx", domain=self.domain, metadata=metadata)

        # Functional Spaces for the wave function and pressure
        self.P2 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), degree_psi)
        self.mixEl = ufl.MixedElement(self.P2, self.P2)
        self.V = FunctionSpace(self.domain, self.mixEl)

        self.P1 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), 1)
        self.Q  = FunctionSpace(self.domain, self.P1)

        # Definition of the trial and test functions for the wave functions
        (self.psi1, self.psi2) = ufl.TrialFunctions(self.V)
        (self.v1,   self.v2)   = ufl.TestFunctions(self.V)
        
        # Initialisation
        # psiTilde = Function(self.V).copy()
        (self.psiTilde_1, self.psiTilde_2) = (Function(self.V.sub(0).collapse()[0]).copy(), Function(self.V.sub(1).collapse()[0]).copy())

        self.pressure = Function(self.Q).copy()

    def assemble(self, hbar: float, direct : bool = False):
        """
        
        This function assembles the forms of the weak formulations and creates the linear structures for the solution (matrix, rhs and solver).

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
        self.a = form( inner(self.psi1, self.v1) * self.dx +
                       inner(self.psi2, self.v2) * self.dx )
        
        # Defining rhs form
        self.L = form( inner(self.psiTilde_1 * ufl.exp( - 1.0 * 1j * self.pressure / self.hbar), self.v1) * self.dx +
                       inner(self.psiTilde_2 * ufl.exp( - 1.0 * 1j * self.pressure / self.hbar), self.v2) * self.dx )
        
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
        self.solution = Function(self.V).copy()
        
    def solve(self, psiTilde_1: dolfinx.fem.Function, psiTilde_2: dolfinx.fem.Function, pressure: dolfinx.fem.Function):
        """
        
        This function computes the shifted wave function.
        
        Parameters
        ----------
        psiTilde_1 : dolfinx.fem.Function
            First component of the wave function :math:`\\tilde{\psi}_1`.
        psiTilde_2 : dolfinx.fem.Function
            Second component of the wave function :math:`\\tilde{\psi}_2`.
        pressure : dolfinx.fem.Function
            Pressure computed by the projection step.
            
        """
        
        # Assign the wave function and pressure to the structures of the class
        with psiTilde_1.vector.localForm() as loc_, self.psiTilde_1.vector.localForm() as loc_n:
            loc_.copy(loc_n)
        with psiTilde_2.vector.localForm() as loc_, self.psiTilde_2.vector.localForm() as loc_n:
            loc_.copy(loc_n)
        with pressure.vector.localForm() as loc_, self.pressure.vector.localForm() as loc_n:
            loc_.copy(loc_n)
        
        # Assembling the rhs vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)

        # Solving the linear system
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()

        return (self.solution.sub(0).collapse(), self.solution.sub(1).collapse())