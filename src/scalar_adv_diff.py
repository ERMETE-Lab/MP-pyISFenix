import numpy as np
import dolfinx
from dolfinx import fem
from dolfinx.fem import Function, FunctionSpace, dirichletbc, locate_dofs_topological, form, assemble_scalar
import ufl
from ufl import grad, div, nabla_grad, inner, dot, real, FacetNormal
from petsc4py import PETSc


class adv_diffusion():
    """
    
    This class implements the solution of the advection-diffusion of a passive scalar in a flow field:
    
    .. math::
        \\frac{\partial T}{\partial t } + \mathbf{u}\cdot \\nabla T-\\alpha\Delta T=0
        
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
        
        # Assign markers - Dirichlet, Neumann, Robin
        self.bound_names = dict()
        self.bound_marks = dict()
        
        for type_ in ['dirichlet', 'neumann', 'robin']:
            if boundary_marks[type_] is not None:
                self.bound_names[type_] = list(boundary_marks[type_].keys())
                self.bound_marks[type_] = [boundary_marks[type_][name] for name in self.bound_names[type_]]
            else:
                self.bound_names[type_] = None
                self.bound_marks[type_] = None

        # Functional Spaces for the temperature and the velocity
        self.P1 = ufl.FiniteElement("Lagrange", self.domain.ufl_cell(), 1)
        self.Q  = FunctionSpace(self.domain, self.P1)

        self.P1vec = ufl.VectorElement("Lagrange", self.domain.ufl_cell(), 1)
        self.W = FunctionSpace(self.domain, self.P1vec)

        # Definition of the trial and test functions for the velocity
        self.T = ufl.TrialFunction(self.Q)
        self.phi = ufl.TestFunction(self.Q)
        
        # Initialisation
        self.Told = Function(self.Q).copy()
        self.uTilde = Function(self.W).copy()
        
    def assemble(self, bc_values: dict, dt = 1e-2, alpha = 1e-3, T_0=0., direct : bool = False):
        """
        
        
        This function assigns the boundary conditions, assembles the forms of the weak formulations and creates the linear structures for the solution (matrix, rhs and solver).

        Parameters
        ----------
        bc_values: dictionary
            Dictionary with the boundary conditions to assign.
        dt : float (Default: 1e-2)
            Time step.
        alpha : float (Default: 1e-3)
            Thermal Diffusivity.
        T_0 : float (Default: 0.)
            Initial value of the temperature.
        direct : boolean (Default: False)
            Boolean variable to identify if a direct solver must be used.
        """
        
        # Defining the diffusivity
        self.alpha = fem.Constant(self.domain, PETSc.ScalarType(alpha))

        # Assign time step size
        self.dt = PETSc.ScalarType(dt)
        
        # Initialise the temperature to T_0
        self.Told.interpolate(lambda x: T_0 + 0.0 * x[0])
        # self.Told.x.set(T_0)
        
        # Defining lhs form
        self.lhs = (inner(self.T / self.dt, self.phi) +
                    inner(dot(self.uTilde, grad(self.T)), self.phi) +
                    inner(self.alpha * grad(self.T), grad(self.phi)) 
                    ) * self.dx
        
        # Defining rhs form
        self.rhs = inner(self.Told / self.dt, self.phi) * self.dx

        self.bcs = []
        # Assigning boundary conditions - Dirichlet
        if self.bound_names['dirichlet'] is not None:
            for ii, name in enumerate(self.bound_names['dirichlet']):
                self.bcs.append(dirichletbc(bc_values['dirichlet'][name], locate_dofs_topological(self.Q, self.fdim, self.ft.find(self.bound_marks['dirichlet'][ii])), self.Q))

        # Assigning boundary conditions - Neumann
        if self.bound_names['neumann'] is not None:
            for ii, name in enumerate(self.bound_names['neumann']):
                self.rhs -= inner(bc_values['neumann'][name], self.phi) * self.ds(self.bound_marks['neumann'][ii])
                
        # Assigning boundary conditions - Robin: 0-th element is the HTC, 1-st is the T_infty
        if self.bound_names['robin'] is not None:
            for ii, name in enumerate(self.bound_names['robin']):
                self.lhs += inner(bc_values['robin'][name][0] * self.T, self.phi) * self.ds(self.bound_marks['robin'][ii])
                self.rhs += inner(bc_values['robin'][name][0] * bc_values['robin'][name][1], self.phi) * self.ds(self.bound_marks['robin'][ii])

        # Defining variational forms
        self.a = form(self.lhs)
        self.L = form(self.rhs)
        
        # Create matrix and vector
        self.A = fem.petsc.create_matrix(self.a)
        self.b = fem.petsc.create_vector(self.L)

        # Define solver 
        self.solver = PETSc.KSP().create(self.domain.comm)
        self.solver.setOperators(self.A)
        if direct:
            self.solver.setType(PETSc.KSP.Type.PREONLY)
            self.solver.getPC().setType(PETSc.PC.Type.LU)
        else:
            self.solver.setType(PETSc.KSP.Type.BCGS)
            # self.solver.getPC().setType(PETSc.PC.Type.ILU)
            self.solver.getPC().setType(PETSc.PC.Type.JACOBI)
        
        # Define solution functions
        self.solution = Function(self.Q).copy()
        
    def solve(self, uTilde: dolfinx.fem.Function):
        """
        
        This function solves the advection-diffusion problem.
        
        Parameters
        ----------
        uTilde : dolfinx.fem.Function
            Velocity field coming from the ISF step.
            
        """
        
        # Assign the velocity to the structures of the class
        with uTilde.vector.localForm() as loc_, self.uTilde.vector.localForm() as loc_n:
            loc_.copy(loc_n)

        # Assembling the matrix
        self.A.zeroEntries()
        fem.petsc.assemble_matrix(self.A, self.a, self.bcs)
        self.A.assemble()

        # Assembling the vector
        with self.b.localForm() as loc:
            loc.set(0)
        fem.petsc.assemble_vector(self.b, self.L)
        fem.petsc.apply_lifting(self.b, [self.a], [self.bcs])
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        fem.petsc.set_bc(self.b, self.bcs)

        # Solving the linear system
        self.solver.solve(self.b, self.solution.vector)
        self.solution.x.scatter_forward()

        with self.solution.vector.localForm() as loc_, self.Told.vector.localForm() as loc_n:
            loc_.copy(loc_n)

        return self.solution