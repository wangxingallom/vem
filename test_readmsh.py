import numpy as np
import sys
import meshio
import matplotlib.pyplot as plt
from fealpy.mesh import TriangleMesh 


point, cell, _, _, _ = meshio.read(sys.argv[1])

node = point[:, 0:2]
cell = cell['triangle']

isUsingNode = np.zeros(node.shape[0], dtype=np.bool)
isUsingNode[cell] = True

NN = isUsingNode.sum()
idxmap = np.zeros(node.shape[0], dtype=np.int32)
idxmap[isUsingNode] = range(NN)
cell = idxmap[cell]
node = node[isUsingNode] 

tmesh = TriangleMesh(node, cell)
NN = tmesh.number_of_nodes()
NC = tmesh.number_of_cells()

fig = plt.figure()
axes = fig.gca()
tmesh.add_plot(axes)
tmesh.find_node(axes)
plt.show()
