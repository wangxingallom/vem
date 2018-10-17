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



pde = CosCosData()
tmesh = unitcircledomainmesh(0.4, meshtype='tri')
pmesh = tri_to_polygonmesh(tmesh,2) 

maxit = 3
Ndof = np.zeros((maxit,), dtype=np.int)
ps = [1]
q  = [5]
#l2
integrator = tmesh.integrator(5)
fem = PoissonFEMModel(pde, tmesh, 1, integrator)
vem = PoissonVEMModel(pde, pmesh, 1, integrator)
fem.solve()
vem.solve()
uv = vem.uh
uf = fem.uh
uI = fem.uI
u = vem.uI
isBdNode = tmesh.ds.boundary_node_flag()
isBdNode1 = pmesh.ds.boundary_node_flag()
e_f = uf[~isBdNode] - uI[~isBdNode]
e_v = uv[~isBdNode1] - u[~isBdNode1]
#L2
femspace = LagrangeFiniteElementSpace(tmesh, 1) 
uhv = femspace.function()
NN = tmesh.number_of_nodes()
uhv[:] = vem.uh[:NN]
u = fem.pde.solution
L2 = fem.integralalg.L2_error(u, uhv)
#H1
gu = fem.pde.gradient
guh = uhv.grad_value
H1 = fem.integralalg.L2_error(gu, guh)
#error
errorMatrix = np.zeros((2, maxit), dtype=np.float)
j=0
i=0    	
errorMatrix[j, i] = np.sqrt(np.mean(e_f**2))
errorMatrix[j+1, i] = np.sqrt(np.mean(e_v**2))
errorMatrix[j, i+1] = fem.get_L2_error()
errorMatrix[j+1, i+1] = L2
errorMatrix[j, i+2] = fem.get_H1_error()
errorMatrix[j+1, i+2] = H1

print(errorMatrix)



