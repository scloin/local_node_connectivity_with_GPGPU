#include <stdio.h>
#include <stdlib.h>
//#include <random>
#include <list>
#include <cstring>

struct pool{
    int source; int target; int *d_label; int *d_visited;
    int* d_edges; int* d_dest; int* h_label; int* h_visited; 
    int numVertex; int* h_returned; int *d_frontier;
    int *d_c_frontier_tail; int *d_p_frontier_tail;
};

void insert_frontier(int source, int* frontier, int* frontier_tail);

int degree(int* dest,int* edges,int source);

void _bidirectional_pred_succ(int* dest,int* edges,int source,int target,int* label_S,int* label_T,int* returned, int ENUM);

int path(int* exclude, int num, int* h_dest, int* h_edges, int state, int ENUM);

void compute_test(int* h_dest,int * h_edges, pool P);