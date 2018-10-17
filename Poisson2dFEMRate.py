import sys

import numpy as np  
import matplotlib.pyplot as plt

from fealpy.pde.poisson_model_2d import CosCosData
from fealpy.fem.PoissonFEMModel import PoissonFEMModel
from fealpy.vem.PoissonVEMModel import PoissonVEMModel
from fealpy.tools.show import showmultirate
from fealpy.mesh.simple_mesh_generator import unitcircledomainmesh, tri_to_polygonmesh
from mpl_toolkits.mplot3d import Axes3D
from fealpy.functionspace.lagrange_fem_space import LagrangeFiniteElementSpace



maxit = 3
pde = CosCosData()
h = 0.1
n = 4
tmesh = unitcircledomainmesh(h, meshtype='tri')
pmesh = tri_to_polygonmesh(tmesh, n) 
errorMatrix = np.zeros((8, maxit), dtype=np.float)
Ndof = np.zeros((maxit,), dtype=np.int)
for i in range(maxit):
    fem = PoissonFEMModel(pde, tmesh, 1, tmesh.integrator(5))
    vem = PoissonVEMModel(pde, pmesh, 1, pmesh.integrator(5))
    fem.solve()
    vem.solve()
    u0 = vem.uh

    u1 = fem.uh
    uI = fem.uI
    u2 = fem.femspace.function()
    N = u2.shape[0]
    u2[:] = u0[:N]

    isBdNode = tmesh.ds.boundary_node_flag()
    errorMatrix[0, i] = np.sqrt(np.mean((uI - u1)**2))
    errorMatrix[1, i] = np.sqrt(np.mean((uI - u2)**2))

    errorMatrix[2, i] = fem.get_L2_error()
    errorMatrix[3, i] = fem.get_L2_error(u2)

    errorMatrix[4, i] = fem.get_H1_error()
    errorMatrix[5, i] = fem.get_H1_error(u2)

    errorMatrix[6, i] = np.max(np.abs(uI - u1)) 
    errorMatrix[7, i] = np.max(np.abs(uI - u2)) 
    
    if i < maxit -1:
        h /= 2
        tmesh = unitcircledomainmesh(h, meshtype='tri')
        pmesh = tri_to_polygonmesh(tmesh, n) 

print(errorMatrix)

