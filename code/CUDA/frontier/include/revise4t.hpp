#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <list>
#include <common/common.h>
#include <common/common_string.h> 
#include <thread>

struct pool{
    int source; int target; 
    int *d_label; int *d_visited; int *d_label1; int *d_visited1;
    int* d_edges; int* d_dest; 
    int* h_label; int* h_visited;  int *h_label1; int *h_visited1;
    int numVertex; int* h_returned; 
    int *d_frontier;
    int *d_c_frontier_tail; int *d_p_frontier_tail;
    cudaStream_t stream;
};

void BFS_host(pool P); 

__global__ void BFS_Bqueue(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited, int* label1, int* visited1); 

__global__ void memset_kernel(int S, int T, int* d_label,int* d_frontier,int numVertex,int* d_p_frontier_tail); 

int degree(int* dest,int* edges,int source); 

int path(int* exclude, int num, int* h_dest, int* h_edges, int state, int elen); 

pool init_pool(int elen, int dlen, int* devmem, int* d_edges);

void compute(int* h_dest,int * h_edges, pool P,FILE* fp1);

__global__ void BFS_Bqueue(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited);

/////////////////////////////

void BFS_host2(pool P);

__global__ void memset_kernel2(int* d_label, int* d_visited, int S, int T, int NUM);

__global__ void BFS_less_atomic(int* d_edges, int* d_dest, int* d_label, int* d_visited, int numVertex, int level);

////////////////////////////

void compute_test(int* h_dest,int * h_edges, pool P);