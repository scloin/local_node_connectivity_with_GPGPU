
import time
import networkx as nx
import sys
from networkx.algorithms.approximation import local_node_connectivity
from networkx.algorithms.connectivity.connectivity import local_node_connectivity as lc

f1 = open("./result/pyver.txt", 'w')
#f2 = open("./result/test0.txt", 'w')
f = open("./graph/graph.txt", 'r')

G=nx.Graph()

destN=int(f.readline())
edgeN=int(f.readline())
dest=f.readline().split(' ')
dest=list(map(int,dest[0:-1]))
edge=f.readline().split(' ')
edge=list(map(int,edge))
for i in range(edgeN-1):
    for j in range(edge[i+1]-edge[i]):
        G.add_edge(i,dest[edge[i]+j])
N=edgeN-1
for i in range(N):
    for j in range(i+1,N):
        #print(i)
        #print(j)
        #print(edge[i])
        #print(edge[j])
        if (type(nx.neighbors(G,i))==int)or(type(nx.neighbors(G,j))==int):
            f1.write("[%d, %d] %d\n" %(i,j,0))
        else:
            f1.write("[%d, %d] %d\n" %(i,j,local_node_connectivity(G,i,j)))
        #f2.write("[%d, %d] %d\n" %(i,j,lc(G,i,j)))
#f2.close()
f1.close()
f.close()
