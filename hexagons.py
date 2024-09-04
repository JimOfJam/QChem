import networkx as nx
import numpy as np
import matplotlib.pyplot as plt 
import scipy
import timeit 

def label(n):
    i = int(np.floor(np.sqrt(n/6)))
    j = int(np.floor((n-6*i**2)/(1+2*i)))
    k = int((n-6*i**2) % (1+2*i))
    return (i,j,k)

def labinv(triple):
    (i,j,k) = triple
    return 6*i**2 + j*(1+2*i) + k

# Construct individual rings
rings = 20
nodes = range(6*(rings)**2)
edges = []
for i in range(rings):
    for j in range(6):
        for k in range(2*i):
            edges.append((labinv((i,j,k)),labinv((i,j,k+1))))
        edges.append((labinv((i,j,2*i)),labinv((i,(j+1) % 6,0))))

# Connect rings together
for i in range(rings-1):
    for j in range(6):
        for k in range(1+2*i):
            if k % 2 == 0:
                edges.append((labinv((i,j,k)),labinv((i+1,j,k+1))))


G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)
nx.draw_kamada_kawai(G,node_size=5)
plt.show()

print("Carbon atoms: " + str(6*rings**2))
alpha = -0.6
beta = -0.1
adjacency = nx.adjacency_matrix(G)
F = alpha*scipy.sparse.identity(6*(rings)**2) + beta*adjacency

# Calculate the eigenvalues of F
start = timeit.default_timer()
spectrum = np.linalg.eigvals(F.todense())
end = timeit.default_timer()
print("Time taken to compute eigenvalues: " + str(end-start) + " seconds.")

plt.hist(spectrum,range=(min(spectrum),max(spectrum)),bins=250)
plt.show()