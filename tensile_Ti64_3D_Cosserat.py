#
# .. raw:: html
#
#  <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><p align="center"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a></p>
#
# .. _PlasticityMFront:
#
# ======================================================================
# Elasto-plastic analysis implemented using the `MFront` code generator
# ======================================================================
#
# This numerical tour has been written in collaboration with **Thomas Helfer** (thomas.helfer@cea.fr), MFront's main developper.
#
# -------------
# Introduction
# -------------
#
# This example is concerned with the incremental analysis of an elasto-plastic
# von Mises material. The behaviour is implemented using the code generator
# tool `MFront` [HEL2015]_ (see http://tfel.sourceforge.net). The behaviour integration 
# is handled by the `MFrontGenericInterfaceSupport` project (see
# https://github.com/thelfer/MFrontGenericInterfaceSupport) which
# allows to load behaviours generated by `MFront` and handle the behaviour integration.
# Other information concerning FEniCS binding in `MFront` can be found at
# https://thelfer.github.io/mgis/web/FEniCSBindings.html, in particular binding
# through the C++ interface and inspired by the `Fenics Solid Mechanics project <https://bitbucket.org/fenics-apps/fenics-solid-mechanics>`_
# is discussed.
#
# The considered example is exactly the same as the pure FEniCS tutorial :ref:`vonMisesPlasticity` so
# that both implementations can be compared. Many implementation steps are
# common to both demos and will not be discussed again here, the reader can refer
# to the previous tutorial for more details.  The sources can be downloaded from 
# :download:`plasticity_mfront.zip`.
#
# Let us point out that a pure FEniCS implementation was possible **only for this specific constitutive law**,
# namely von Mises plasticity with isotropic linear hardening, since the return mapping step
# is analytical in this case and can be expressed using simple UFL operators. The use of `MFront`
# can therefore make it possible to consider more complex material laws. Indeed, one only has to change 
# the names of the behaviour and library generated by `MFront` to change the material behaviours, which 
# can be arbitraly complex (although limited to small strains). This will be illustrated in forthcoming tutorials. 
# Interested users may have a look at the examples of the `MFront` gallery to have a small overview of 
# `MFront` abilities: http://tfel.sourceforge.net/gallery.html
#
# -------------
# Prerequisites
# -------------
# In order to run this numerical tour, you must first install the TFEL library on which `MFront`
# is built as well as MGIS (`MFrontGenericInterfaceSupport`). Please note that development
# versions are used for now.
#
# To proceed with the installation, we provide the following installation shell script :download:`install.sh`. 
# Before running it, please install the necessary prerequisites mentioned in the script, e.g. on Ubuntu run::

# > sudo apt-get install cmake libboost-all-dev g++ gfortran
# > sudo apt-get install git libqt5svg5-dev qtwebengine5-dev
# > sudo apt-get install python3-matplotlib
#
# Once the installation script finished installing TFEL and MGIS, you can use MGIS by running the
# following command::

# > source <PREFIX>/codes/mgis/master/install/env.sh
#
# in which ``<PREFIX>`` is the installation directory you have chosen.
#
# As recalled later in the tutorial, the MFront behaviour file must first be compiled before running 
# this Python script. The compilation command is::

# > mfront --obuild --interface=generic IsotropicLinearHardeningPlasticity.mfront
#
# -----------------
# Problem position
# -----------------
#
# The present implementation will heavily rely on ``Quadrature`` elements to represent
# previous and current stress states, internal state variables (cumulated
# equivalent plastic strain :math:`p` in the present case) as well as components
# of the tangent stiffness matrix. The use of ``quadrature`` representation is therefore
# still needed to circumvent the known issue with ``uflacs``::

from __future__ import print_function
from dolfin import *
import numpy as np
parameters["form_compiler"]["representation"] = 'quadrature'
import warnings
from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning
warnings.simplefilter("once", QuadratureRepresentationDeprecationWarning)
from ufl import Identity, PermutationSymbol, indices
import logging
logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)

# In order to bind `MFront` behaviours with FEniCS Python interface we will first load
# the `MFrontGenericInterfaceSupport` Python package named ``mgis`` and more precisely
# the ``behaviour`` subpackage::

import mgis.behaviour as mgis_bv

# As before, we load a ``Gmsh`` generated mesh representing a portion of a thick cylinder::

#Re, Ri = 1.3, 1.   # external/internal radius
#mesh = Mesh("tensile_dog_bone_specimen.xml")
#facets = MeshFunction("int", mesh, "./tensile_dog_bone_specimen_gmsh:geometrical.xml")

mesh = UnitCubeMesh(1,1,1)
#mesh = Mesh()
'''XDMFFile("mesh.xdmf").read(mesh)
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()-1)
XDMFFile("mf.xdmf").read(mvc, "name_to_read")
facets = cpp.mesh.MeshFunctionSizet(mesh, mvc)'''

#ds = Measure('ds')[facets]
eps_ds = 1E-14

#def boundary_L(x, on_boundary):
#    return on_boundary and near(x[0], -80.0, eps_ds)

def boundaryX0(x, on_boundary):
    return on_boundary and near(x[0], 0.0, eps_ds)

#def boundary_R(x, on_boundary):
#    return on_boundary and near(x[0], 80.0, eps_ds)

def boundaryY0(x, on_boundary):
    return on_boundary and near(x[1], 0.0, eps_ds)

def boundaryZ0(x, on_boundary):
    return on_boundary and near(x[2], 0.0, eps_ds)


# We now define appropriate function spaces. Standard CG-space of degree 2 will still be used 
# for the displacement whereas various Quadrature spaces are considered for:
#
# * stress/strain-like variables (represented here as 4-dimensional vector since the :math:`zz` component must be considered in the plane strain plastic behaviour)
# * scalar variables for the cumulated plastic strain
# * the consistent tangent matrix represented here as a tensor of shape 4x4
#
# As in the previous tutorial a ``degree=2`` quadrature rule (i.e. 3 Gauss points)
# will be used. In the end, the total number of Gauss points in the mesh is retrieved
# as it will be required to instantiate `MFront` objects (note that it can be obtained
# from the dimension of the Quadrature function spaces or computed from the number of
# mesh cells and the chosen quadrature degree)::

deg_u = 2
deg_stress = 2
stress_strain_dim = 9
stress_dim = 3

VE = VectorElement("CG", mesh.ufl_cell(), 2)
THETAE = VectorElement("CG", mesh.ufl_cell(), 2)
VT = FunctionSpace(mesh, MixedElement((VE, THETAE)))

# Quadrature space for sigma
We = VectorElement("Quadrature", mesh.ufl_cell(), degree=deg_stress,
                    dim=stress_strain_dim,
		    quad_scheme='default')
W = FunctionSpace(mesh, We)
# Quadrature space for p
W0e = FiniteElement("Quadrature", mesh.ufl_cell(), degree=deg_stress, quad_scheme='default')
W0 = FunctionSpace(mesh, W0e)
# Quadrature space for tangent matrix
Wce = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg_stress,
                     shape=(stress_strain_dim, stress_strain_dim),
                     quad_scheme='default')
WC = FunctionSpace(mesh, Wce)

# get total number of gauss points
ngauss = W0.dim()

#print("# GAUSS points: ", ngauss)

# Various functions are defined to keep track of the current internal state (stresses, 
# current strain estimate, cumulative plastic strain and cpnsistent tangent matrix)
# as well as the previous displacement, current displacement estimate and current iteration correction::

# Define functions based on Quadrature spaces
sig = Function(W, name="Current stresses")
sig_old = Function(W)
couple_sig = Function(W, name="Current couple stresses")
couple_sig_old = Function(W)
Eps1 = Function(W, name="Current strain increment estimate at the end of the end step")
Kappa1 = Function(W, name="Current wryness increment estimate at the end of the end step")
Ct = Function(WC, name="Stress-Strain consistent tangent operator")
Dt = Function(WC, name="CupleStress-Wryness consistent tangent operator")
p = Function(W0, name="Cumulative plastic strain")

vt = Function(VT)
(u, theta) = vt.split()

vt1 = Function(VT)
(u1, theta1) = vt1.split()

dvt = Function(VT)
(du, dtheta) = dvt.split()

(u_v, t_v) = TrialFunctions(VT)

(u_, t_) = TestFunctions(VT)

'''u  = Function(V, name="Displacement at the beginning of the time step")
u1 = Function(V, name="Current displacement estimate at the end of the end step")
du = Function(V, name="Iteration correction")
v = TrialFunction(V)
u_ = TestFunction(V)'''


#print("P get local : ", dir(p.get_local()))
#print("P get local : ", p.get_local())

# ----------------------------------------------------
# Material constitutive law definition using `MFront`
# ----------------------------------------------------
#
# We now define the material. First let us make a rapid description of some classes 
# and functions introduced by the `MFrontGenericInterface` project that will be helpful for this
# tutorial:
#
# * the ``Behaviour`` class handles all the information about a specific
#   `MFront` behaviour. It is created by the ``load`` function which takes
#   the path to a library, the name of a behaviour and a modelling
#   hypothesis.
# * the ``MaterialDataManager`` class handles a bunch of integration points.
#   It is instantiated using an instance of the ``Behaviour`` class and the
#   number of integration points [#]_. The ``MaterialDataManager`` contains
#   various interesting members:
#  
#   - ``s0``: data structure of the ``MaterialStateManager`` type which stands for the material state at the beginning of the time step.
#   - ``s1``: data structure of the ``MaterialStateManager`` type which stands for the material state at the end of the time step.
#   - ``K``: a ``numpy`` array containing the consistent tangent operator at each integration points.
#      
# * the ``MaterialStateManager`` class describe the state of a material. The
#   following members will be useful in the following:
#    
#     - ``gradients``: a numpy array containing the value of the gradients
#       at each integration points. The number of components of the
#       gradients at each integration points is given by the
#       ``gradients_stride`` member.
#     - ``thermodynamic_forces``: a numpy array containing the value of the
#       thermodynamic forces at each integration points. The number of
#       components of the thermodynamic forces at each integration points
#       is given by the ``thermodynamic_forces_stride`` member.
#     - ``internal_state_variables``: a numpy array containing the value of the
#       internal state variables at each integration points. The number of
#       internal state variables at each integration points is given by the
#       ``internal_state_variables_stride`` member.
#      
# * the ``setMaterialProperty`` and ``setExternalStateVariable`` functions can
#   be used to set the value a material property or a state variable
#   respectively.
# * the ``update`` function updates an instance of the
#   ``MaterialStateManager`` by copying the state ``s1`` at the end of the
#   time step in the state ``s0`` at the beginning of the time step.
# * the ``revert`` function reverts an instance of the
#   ``MaterialStateManager`` by copying the state ``s0`` at the beginning of
#   the time step in the state ``s1`` at the end of the time step.
# * the ``integrate`` function triggers the behaviour integration at each
#   integration points. Various overloads of this function exist. We will
#   use a version taking as argument a ``MaterialStateManager``, the time
#   increment and a range of integration points.
#
# In the present case, we compute a plane strain von Mises plasticity with isotropic
# linear hardening. The material behaviour is implemented in the :download:`IsotropicLinearHardeningPlasticity.mfront` file
# which must first be compiled to generate the appropriate librairies as follows::

# > mfront --obuild --interface=generic IsotropicLinearHardeningPlasticity.mfront
#
# We can then setup the ``MaterialDataManager``::

# Defining the modelling hypothesis
h = mgis_bv.Hypothesis.TRIDIMENSIONAL
# Loading the behaviour        
b = mgis_bv.load('src/libBehaviour.so','CosseratIsotropicLinearHardeningPlasticity',h)
# Setting the material data manager
m = mgis_bv.MaterialDataManager(b, ngauss)
# elastic parameters for Cosserat strain
E = 107e3
llambda = 107e3 # In [MPa]
mu = 70e3
mu_c = 100e3
# elastic parameters for Cosserat wryness
alpha = 0. # In [MPa]
beta  = 30e3
gamma = 30e3
# yield strength
sig0 = 800.
Et = E/100.
# hardening slope
H = E*Et/(E-Et)
# Cosserat plasticity coupling modules
a_1 = 1.0
a_2 = 1.0
b_1 = 1.0
b_2 = 1.0

for s in [m.s0, m.s1]:
    mgis_bv.setMaterialProperty(s, "lambda", llambda)
    mgis_bv.setMaterialProperty(s, "mu", mu)
    mgis_bv.setMaterialProperty(s, "mu_c", mu_c)
    mgis_bv.setMaterialProperty(s, "alpha", alpha)
    mgis_bv.setMaterialProperty(s, "beta", beta)
    mgis_bv.setMaterialProperty(s, "gamma", gamma)
    mgis_bv.setMaterialProperty(s, "a_1", a_1)
    mgis_bv.setMaterialProperty(s, "a_2", a_2)
    mgis_bv.setMaterialProperty(s, "b_1", b_1)
    mgis_bv.setMaterialProperty(s, "b_2", b_2)
    mgis_bv.setMaterialProperty(s, "HardeningSlope", H)
    mgis_bv.setMaterialProperty(s, "YieldStrength", sig0)
    mgis_bv.setExternalStateVariable(s, "Temperature", 293.15)

# Boundary conditions and external loading are defined as before along with the 
# analytical limit load solution::

bc = [DirichletBC(VT.sub(0), Constant((0,0,0)), boundaryX0), DirichletBC(VT.sub(1), Constant((0,0,0)), boundaryX0)] #Clamped Dispalcements and Rotations on the surface at x=0
bc = [DirichletBC(VT.sub(0).sub(1), Constant(0), boundaryY0)] #Fixed Dispalcement along Y on the surface at y=0
bc = [DirichletBC(VT.sub(0).sub(2), Constant(0), boundaryZ0)] #Fixed Dispalcement along Z on the surface at z=0
#bc = [DirichletBC(V.sub(0), 0, facets, 38), DirichletBC(V.sub(1), 0, facets, 38),  DirichletBC(V.sub(2), 0, facets, 38)]
#bc = [DirichletBC(V, as_vector((0.0,0.0,0.0)), facets, 1)]
eps_load = 1e-5
#loadside = Expression("x[0] < 80.0 + eps && x[0] > 80.0 - eps ? 1. : 0.", eps=eps_load, degree=2)
loadside = Expression("x[0] < 1.0 + eps && x[0] > 1.0 - eps ? 1. : 0.", eps=eps_load, degree=2)

n = FacetNormal(mesh)
#q_lim = float(2/sqrt(3)*ln(Re/Ri)*sig0)
q_lim = 1050 #[MPa]
loading = Expression("q*t", q=q_lim, t=0, degree=2)

#print(dir(ds), ds(2000).subdomain_id(), ds(2000).subdomain_data().size(), ds(2000).integral_type())

def F_ext(v):
    return loadside * loading * n[0] * v[0] * ds

# --------------------------------------------
# Global problem and Newton-Raphson procedure
# --------------------------------------------
#
# Before writing the global Newton system, we first define the strain measure 
# function ``eps_MFront`` consistent with the `MFront` conventions (see 
# http://tfel.sourceforge.net/tensors.html for details). We also define the custom
# quadrature measure and the projection function onto Quadrature spaces::

def eps_MFront(v_u, v_theta):
    e = grad(v_u) +  as_matrix(((0,v_theta[2],-v_theta[1]),(-v_theta[2],0,v_theta[0]),(v_theta[1],-v_theta[0],0)))
    return as_vector([e[0, 0], e[1, 1], e[2,2], e[0,1], e[1,0], e[0,2], e[2,0], e[1,2], e[2,1]])

def kappa_MFront(v_theta):
    k = grad(v_theta)
    return as_vector([k[0, 0], k[1, 1], k[2,2], k[0,1], k[1,0], k[0,2], k[2,0], k[1,2], k[2,1]])


metadata = {"quadrature_degree": deg_stress, "quadrature_scheme": "default"}
dxm = dx(metadata=metadata)

def local_project(v, V, u=None):
    """ 
    projects v on V with custom quadrature scheme dedicated to
    FunctionSpaces V of `Quadrature` type
        
    if u is provided, result is appended to u
    """
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dxm
    b_proj = inner(v, v_)*dxm
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

def exact_project(v, V, u):
    """ 
    projects v on V with custom quadrature scheme dedicated to
    FunctionSpaces V of `Quadrature` type
        
    if u is provided, result is appended to u
    """
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dxm
    b_proj = inner(v, v_)*dxm
    solve(a_proj == b_proj, u)



# The bilinear form of the global problem is obtained using the consistent tangent
# matrix ``Ct`` and the `MFront` strain measure, whereas the right-hand side consists of
# the residual between the internal forces associated with the current
# stress state ``sig`` and the external force vector. ::

a_Newton = inner(eps_MFront(u_v, t_v), dot(Ct, eps_MFront(u_, t_)))*dxm + inner(kappa_MFront(t_v), dot(Dt, kappa_MFront(t_)))*dxm
res = -inner(eps_MFront(u_, t_), sig)*dxm - inner(kappa_MFront(t_), couple_sig)*dxm + F_ext(u_)

# Before defining the Newton-Raphson loop, we set up the output file and appropriate
# FunctionSpace (here piecewise constant) and Function for output of the equivalent
# plastic strain since XDMF output does not handle Quadrature elements::

file_results = XDMFFile("plasticity_results.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
P0 = FunctionSpace(mesh, "DG", 0)
p_avg = Function(P0, name="Plastic strain")

# The tangent stiffness is also initialized with the elasticity matrix:

# one integrates the behaviour over the time step and computes an elastic stiffness
it = mgis_bv.IntegrationType.PredictionWithElasticOperator
#print("Numer of Gauss Points m.n: ", m.n)
#m.n numer of intergation points
mgis_bv.integrate(m, it, 0, 0, m.n);
tangent_operators = m.K.flatten()
#print("m.K_stride:",m.K_stride)
#print("K:",m.K)
Ct.vector().set_local(tangent_operators[0:m.n*stress_strain_dim**2-1])
Ct.vector().apply("insert")
# getting the coupled tangent operator
Dt.vector().set_local(tangent_operators[m.n*stress_strain_dim**2:])
Dt.vector().apply("insert")

# The main difference with respect to the pure FEniCS implementation of the previous
# tutorial is that `MFront` computes the current stress state and stiffness matrix
# (``integrate`` method) based on the value of the total strain which is computed 
# from the total displacement estimate ``u1``. The associated strain is projected 
# onto the appropriate Quadrature function space so that its array of values at all 
# Gauss points is given to `MFront` via the ``m.s1.gradients`` attribute. The flattened
# array of stress and tangent stiffness values are then used to update the current 
# stress and tangent stiffness variables. The cumulated plastic strain is also
# retrieved from the ``internal_state_variables`` attribute (:math:`p` being the last 
# column in the present case). At the end of the iteration loop, the material 
# behaviour and the previous displacement variable are updated::

Nitermax, tol = 1, 1e-5  # parameters of the Newton-Raphson procedure
Nincr = 100
load_steps = np.linspace(0, 1., Nincr+1)[1:]
results = np.zeros((Nincr+1, 2))

# creating a thread pool for parallel integration
#pool = mgis_bv.ThreadPool(4)

for (i, t) in enumerate(load_steps[0:1]):
    loading.t = t
    A, Res = assemble_system(a_Newton, res, bc)
    #print(A, "RES: ", Res.get_local(), " size: ", Res.size(), " value: ", Res.norm("l2"))
    nRes0 = Res.norm("l2")
    nRes = nRes0
    print("nRes:",nRes)
    vt1.assign(vt)
    #theta1.assign(theta)
    print("Increment:", str(i+1), " , t:", t)
    niter = 0
    while nRes > tol and niter < Nitermax:
        print("    Residual ratio: ", nRes/nRes0, ", nRes0: ", nRes0)
        solve(A, dvt.vector(), Res, "mumps")
        #the current estimate of the displacement at the end of the time step
        vt1.assign(vt1+dvt)
        #u.assign(u1+du)
        #theta1.assign(theta1+dtheta)
        # compute the current estimate of the strain at the end of the
        # time step using `MFront` conventions
        (u1, theta1) = vt1.split()
        exact_project(eps_MFront(u1, theta1), W, Eps1)
        #local_project(eps_MFront(u1, theta1), W, Eps1)
        exact_project(kappa_MFront(theta1), W, Kappa1)
        #local_project(kappa_MFront(theta1), W, Kappa1)
        # copy the strain values to `MGIS`
        m.s1.gradients[:, 0:stress_strain_dim] = Eps1.vector().get_local().reshape((m.n, stress_strain_dim))
        m.s1.gradients[:, stress_strain_dim:] = Kappa1.vector().get_local().reshape((m.n, stress_strain_dim))
        #print("Displacement u: ", du.vector().get_local())
        #print("Rotation theta: ", dtheta.vector().get_local())
        #print("Eto: ", m.s1.gradients[:, :stress_strain_dim])
        #print("Kappa: ", m.s1.gradients[:, stress_strain_dim:])
        # integrate the behaviour
        it = mgis_bv.IntegrationType.IntegrationWithConsistentTangentOperator
        mgis_bv.integrate(m, it, 0, 0, m.n);
        # getting the stress and consistent tangent operator back to
        # the FEniCS world.
        stress_states = m.s1.thermodynamic_forces.flatten()
        sig.vector().set_local(stress_states[0:m.n*stress_strain_dim])
        sig.vector().apply("insert")
        #print("Stress state: ",stress_states[0:m.n*stress_strain_dim])
        # getting the couple stress
        couple_sig.vector().set_local(stress_states[m.n*stress_strain_dim:])
        couple_sig.vector().apply("insert")
        #print("Coupled Stress state: ",stress_states[m.n*stress_strain_dim:])
        # getting the tangent operator
        tangent_operators = m.K.flatten()
        Ct.vector().set_local(tangent_operators[0:m.n*stress_strain_dim**2])
        Ct.vector().apply("insert")
        print("Tangent Operator:", tangent_operators[0:m.n*stress_strain_dim**2])
        # getting the coupled tangent operator
        Dt.vector().set_local(tangent_operators[m.n*stress_strain_dim**2:])
        Dt.vector().apply("insert")
        #print("Coupled Tangent Operator:", tangent_operators[m.n*stress_strain_dim**2:])
        # retrieve cumulated plastic strain values
        p.vector().set_local(m.s1.internal_state_variables[:, -1])
        p.vector().apply("insert")
        # solve Newton system
        A, Res = assemble_system(a_Newton, res, bc)
        nRes = Res.norm("l2")
        print("    Residual at the end:", nRes)
        niter += 1
    # update the displacement for the next increment
    vt.assign(vt1)
    # update the material
    # the update function updates an instance of the MaterialStateManager by copying the state s1 at the end of the time step in the state s0 at the beginning of the time step.
    mgis_bv.update(m)
    
    # postprocessing results
    file_results.write(u, t)
    p_avg.assign(project(p, P0))
    file_results.write(p_avg, t)
    results[i+1, :] = (u(0, 0, 0)[0], t)

import matplotlib.pyplot as plt
plt.plot(results[:, 0], results[:, 1], "-o")
plt.xlabel("Displacement of inner boundary")
plt.ylabel(r"Applied pressure $q/q_{lim}$")
plt.show()

# .. note::
#  Note that we defined the cumulative plastic strain variable :math:`p` in FEniCS
#  only for post-processing purposes. In fact, FEniCS deals only with the global equilibrium
#  whereas `MFront` manages the history of internal state variables, so that this variable
#  would not have been needed if we were not interested in post-processing it.
# 
# -------------------------
# Results and future works
# -------------------------
#
# We can verify that the convergence of the Newton-Raphson procedure is extremely similar
# between the `MFront`-based implementation and the pure FEniCS one, the same number of 
# iterations per increment is obtained along with close values of the residual.
#
# **Total computing time** took approximately:
#
# * 5.9s for the present `MFront` implementation against
# * 6.8s for the previous FEniCS-only implementation
#  
# Several points need to be mentioned regarding this implementation efficiency:
#
# * MGIS can handle parallel integration of the constitutive law which has not been used
#   for the present computation
# * the present approach can be improved by letting MGIS reuse the memory already allocated
#   by FEniCS which will reduce information transfer times and memory consumption
# * extension to large strains is a work in progress
# * this FEniCS/MGIS coupling will make it possible, in a near future, to test in a rapid
#   manner generalized constitutive laws (higher-order theories, phase-field) and/or
#   multiphysics couplings
#   
# ------------
#  References
# ------------
#
# .. [HEL2015] Helfer, Thomas, Bruno Michel, Jean-Michel Proix, Maxime
#  Salvo, Jérôme Sercombe, and Michel Casella. 2015. *Introducing the
#  Open-Source Mfront Code Generator: Application to Mechanical
#  Behaviours and Material Knowledge Management Within the PLEIADES
#  Fuel Element Modelling Platform.* Computers & Mathematics with
#  Applications. <https://doi.org/10.1016/j.camwa.2015.06.027>.
#
# .. [#] Note that an instance of MaterialDataManager keeps a reference to the behaviour 
# which has been used for its initialization: the user must ensure that this behaviour 
# outlives the instance of the MaterialDataManager, otherwise memory corruption may occur.
#
