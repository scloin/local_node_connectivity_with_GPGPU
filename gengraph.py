
import time
import networkx as nx
import sys
from networkx.algorithms.approximation import local_node_connectivity
#import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("120, 27(default)")
    N=20
    S = 20
    P=0.4
    #sys.exit()
else:
    N =int(sys.argv[1])
    P= float(sys.argv[2])
    S= int(sys.argv[3])
    #print(N,S)

f = open("./graph/graph.txt", 'w')

G=nx.gnp_random_graph(N,P,S)

dest,edge=[],[0]
t=0
for i in G.nodes:
    temp=[]
    for j in G.neighbors(i):
        temp.append(j)
        t+=1
    dest+=temp
    edge.append(t)
f.write(str(len(dest)))
f.write(" \n")
f.write(str(len(edge)))
f.write(" \n")
f.write(str(' '.join(map(str,dest))))
f.write(" \n")
f.write(str(' '.join(map(str,edge))))
# print(dest)
# print(edge)
f.close()

#nx.draw(G, pos=nx.circular_layout(G), node_color="red")

#plt.show()