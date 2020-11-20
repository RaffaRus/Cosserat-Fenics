# Import FEniCS library
from dolfin import *
import mgis.fenics as mf
import numpy as np
import ufl



# Parameters
E_el=200000.      #Young modulus
nu_el=0.3       #Poisson ratio
eta_1=0.8       #Cosserat coefficient for Gc
eta_3=5./2.     #Cosserat coefficient for Mc
lc=1/10.         #Cosserat length

N = 20            #number of elements - FE discretization

G_el = E_el/(2.*(1+nu_el))
M_el = E_el*(1.-nu_el)/((1.-2.*nu_el)*(1.+nu_el))
Gc_el = eta_1*G_el
Mc_el = 2*eta_3*G_el*lc**2

# Generate mesh
mesh = BoxMesh(Point(0., 0., 0.), Point(20., 5., 20.), N, N, 1)
# Defines a LagrangeFE of degree 2 for the displacements
Ue = VectorElement("Lagrange", mesh.ufl_cell() ,degree=2)
# Defines a Lagrangie FE of degree 1 for the rotations
Te = VectorElement("Lagrange", mesh.ufl_cell(), degree=1)
# Creates a mixed function space
V = FunctionSpace(mesh, MixedElement([Ue, Te]))
# Define test functions (virtual velocities)
v = Function(V)
(u, theta) = split(v)

# Levi-Civita permutation symbol
d = mesh.geometry().dim()
I = np.eye(d)
A = lambda *ind: np.linalg.det(np.array([[I[p, q] for q in range(d)] for p in ind]))
perm = as_tensor([[[A(i, j, k) for k in range(d)] for j in range(d)] for i in range(d)])
print(perm, dir(perm))
eps = grad(u) + dot(perm, theta)
kappa = grad(theta)

# Define boundary conditions
def bottom(x, on_boundary):
    return on_boundary and near(x[1], 0.)
def top(x, on_boundary):
    return on_boundary and near(x[1], 5.)
def whole(x, on_boundary):
    return on_boundary;


# 80 couple traction on top 
# 0.001 grad microrotation 

#Uimp = Expression((0.,)*5 + ("t*0.001",), t=0.2, degree=0)
Uimp = Expression("t*0.001", t=0.2, degree=0)
bc= [DirichletBC(V, Constant((0.,)*6), bottom), DirichletBC(V.sub(0).sub(2), Constant(0.), whole), DirichletBC(V.sub(1).sub(2), Uimp, top), DirichletBC(V.sub(1).sub(0), Constant(0.), whole), DirichletBC(V.sub(1).sub(1), Constant(0.), whole)]

facets = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
AutoSubDomain(top).mark(facets, 1)

file_results = XDMFFile("results/cosserat_plasticity.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True

mu = E_el/2/(1+nu_el)
mat_prop = {
 'lambda': E_el*nu_el/(1+nu_el)/(1-2*nu_el),
 'mu': mu,
 'mu_c': 100000,
 'alpha': 0.,
 'beta': 76923,
 'gamma': 76923,
 'a_1': 1,
 'a_2': 0,
 'b_1': 1,
 'b_2': 0,
 'HardeningSlope': 0.,
 'YieldStrength': 100.}

material = mf.MFrontNonlinearMaterial("./src/libBehaviour.so",
                                          "ExplicitCosseratIsotropicLinearHardeningPlasticity",
                                          material_properties=mat_prop)
problem = mf.MFrontNonlinearProblem(v, material, quadrature_degree=4, bcs=bc)
problem.register_gradient("TotalStrain", eps)
problem.register_gradient("WrynessTensor", kappa)

print(dir(problem))
problem.solver = PETScSNESSolver()
prm = problem.solver.parameters
prm["method"] = "newtonls"
prm["linear_solver"] = "mumps"
tol = 1e-4
prm["solution_tolerance"] = tol
prm["relative_tolerance"] = tol

problem.dt = 1.
for t in np.linspace(0, 1.0, 20)[1:]:
    Uimp.t = t
    problem.solve(v.vector())

    u = v.sub(0, True)
    u.rename("u", "u")
    theta = v.sub(1, True)
    theta.rename("Theta", "Theta")
    p = problem.get_state_variable("EquivalentPlasticStrain", project_on=("CG", 1))
    el = problem.get_state_variable("ElasticStrain", project_on=("CG", 1))
    kel = problem.get_state_variable("ElasticCurvature", project_on=("CG", 1))
    sig = problem.get_flux("Stress", project_on=("CG", 1))
    couple_sig = problem.get_flux("CoupleStress", project_on=("CG", 1))
    file_results.write(u, t)
    file_results.write(theta, t)
    file_results.write(p, t)
    file_results.write(el, t)
    file_results.write(kel, t)
    file_results.write(sig, t)
    file_results.write(couple_sig, t)

