
import time
#import networkx as nx
import sys
#from networkx.algorithms.approximation import local_node_connectivity
#from networkx.algorithms.connectivity.connectivity import local_node_connectivity as lc

INF = float("inf")

class Graph:
    def __init__(self, destN, edgeN, dest, edge):
        self.destN=destN
        self.edgeN=edgeN
        self.dest=dest
        self.edge=edge
        
    def neighbors(self, node):
        return iter(self.dest[self.edge[node]:self.edge[node+1]])
    
    def degree(self, node):
        return len(self.dest[self.edge[node]:self.edge[node+1]])

def local_node_connectivity(G, source, target, cutoff=None):

    if target == source:
        raise ("source and target have to be different nodes.")

    # Maximum possible node independent paths

    possible = min(G.degree(source), G.degree(target))

    K = 0
    if not possible:
        return K

    if cutoff is None:
        cutoff = INF

    exclude = set()
    for i in range(min(possible, cutoff)):
        try:
            path = _bidirectional_shortest_path(G, source, target, exclude)
            exclude.update(set(path))
            K += 1
        except:
            break

    return K

def _bidirectional_shortest_path(G, source, target, exclude):

    # call helper to do the real work
    results = _bidirectional_pred_succ(G, source, target, exclude)
    pred, succ, w = results

    # build path from pred+w+succ
    path = []
    # from source to w
    while w is not None:
        path.append(w)
        w = pred[w]
    path.reverse()
    # from w to target
    w = succ[path[-1]]
    while w is not None:
        path.append(w)
        w = succ[w]

    return path

def _bidirectional_pred_succ(G, source, target, exclude):
    # does BFS from both source and target and meets in the middle
    # excludes nodes in the container "exclude" from the search
    if source is None or target is None:
        raise (
            "Bidirectional shortest path called without source or target"
        )
    if target == source:
        return ({target: None}, {source: None}, source)

    # handle either directed or undirected

    Gpred = G.neighbors
    Gsucc = G.neighbors

    # predecesssor and successors in search
    pred = {source: None}
    succ = {target: None}

    # initialize fringes, start with forward
    forward_fringe = [source]
    reverse_fringe = [target]

    level = 0

    while forward_fringe and reverse_fringe:
        # Make sure that we iterate one step forward and one step backwards
        # thus source and target will only trigger "found path" when they are
        # adjacent and then they can be safely included in the container 'exclude'
        level += 1
        if not level % 2 == 0:
            this_level = forward_fringe
            forward_fringe = []
            for v in this_level:
                for w in Gsucc(v):
                    if w in exclude:
                        continue
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w] = v
                    if w in succ:
                        return pred, succ, w  # found path
        else:
            this_level = reverse_fringe
            reverse_fringe = []
            for v in this_level:
                for w in Gpred(v):
                    if w in exclude:
                        continue
                    if w not in succ:
                        succ[w] = v
                        reverse_fringe.append(w)
                    if w in pred:
                        return pred, succ, w  # found path

    raise (f"No path between {source} and {target}.")


f1 = open("./result/pyver_test.txt", 'w')
#f2 = open("./result/test0.txt", 'w')
f = open("./graph/graph.txt", 'r')



destN=int(f.readline())
edgeN=int(f.readline())
dest=f.readline().split(' ')
dest=list(map(int,dest[0:-1]))
edge=f.readline().split(' ')
edge=list(map(int,edge))
G=Graph(destN,edgeN,dest,edge)
# G=nx.Graph()
# for i in range(edgeN-1):
#     for j in range(edge[i+1]-edge[i]):
#         G.add_edge(i,dest[edge[i]+j])
N=edgeN-1
for i in range(N):
    for j in range(i+1,N):
        try:
            f1.write("[%d, %d] %d\n" %(i,j,local_node_connectivity(G,i,j)))
            #local_node_connectivity(G,i,j)
        except:
            #f1.write("[%d, %d] %d\n" %(i,j,0))
            pass
        break
    break

        
f1.close()
f.close()