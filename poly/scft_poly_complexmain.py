import sys
import numpy as np
import matplotlib.pyplot as plt
import meshio

from fealpy.functionspace.vem_space import VirtualElementSpace2d 
from SCFTVEMModel import SCFTVEMModel, SCFTParameter
from fealpy.mesh.QuadrangleMesh import QuadrangleMesh
from fealpy.mesh.tree_data_structure import Quadtree
from fealpy.mesh import TriangleMesh, TriangleMeshWithInfinityNode 
from fealpy.mesh.PolygonMesh import PolygonMesh

def complex_mesh(r=5):
    mesh = meshio.read(sys.argv[1])
    node = mesh.points
    node = node[:,0:2]*r
    cell = mesh.cells
    cell = cell['triangle']
    isUsingNode = np.zeros(node.shape[0], dtype=np.bool)
    isUsingNode[cell] = True
    NN = isUsingNode.sum()
    idxmap = np.zeros(node.shape[0], dtype=np.int32)
    idxmap[isUsingNode] = range(NN)
    cell = idxmap[cell]
    node = node[isUsingNode]
    mesh = TriangleMesh(node,cell)
    nmesh = TriangleMeshWithInfinityNode(mesh)
    ppoint, pcell, pcellLocation =  nmesh.to_polygonmesh()
    pmesh = PolygonMesh(ppoint, pcell, pcellLocation)
    NN = pmesh.number_of_nodes()
    NC = pmesh.number_of_cells()
    return pmesh

fieldsType = int(sys.argv[2])
option = SCFTParameter()
option.fA = 0.5
option.chiN = 15

option.maxit = 5000
option.showstep = 2
mesh = complex_mesh(r=20) 
NN = mesh.number_of_nodes()

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.show()

print('The number of mesh:', NN)
# define vemspace
p = 1
vemspace =VirtualElementSpace2d(mesh, p)
# define fields
Ndof = vemspace.number_of_global_dofs()
fields = np.zeros((Ndof, 2))
node = mesh.node
chiN = option.Ndeg*option.chiAB

if fieldsType == 1:
    fields[:, 0] = chiN*(-1 + 2*np.random.rand(1, Ndof))
    fields[:, 1] = chiN*(-1 + 2*np.random.rand(1, Ndof))
elif fieldsType == 2:
    pi = np.pi
    fields[:, 1] = chiN * (np.sin(3*node[:,0]) + np.cos(3*node[:, 1]))
elif fieldsType == 3:
    pi = np.pi
    x = node[:, 0]
    y = node[:, 1]
    def f(k,x, y):
        return np.cos(np.cos(pi/3*k)*x + np.sin(pi/3*k)*y)
    fields[:, 1] = np.cos(f(0,x,y)) + np.cos(f(1,x,y)) + np.cos(f(2,x,y)) + np.cos(f(3,x,y)) + np.cos(f(4,x,y)) + np.cos(f(5,x,y)) 
elif fieldsType == 4:
    x = node[:, 0]
    fields[:, 1] = np.sin(1*x)
    
option.fields = fields
scft = SCFTVEMModel(vemspace, option) 
scft.initialize()
scft.find_saddle_point(datafile='huabandata', file_path='./results')

