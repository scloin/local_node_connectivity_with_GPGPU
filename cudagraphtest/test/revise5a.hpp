#include <stdio.h>
#include <stdlib.h>
#include<iostream>
#include<thread>
#include <cuda_runtime.h>
#include <list>
#include "common.h"
#include "common_string.h"
#include "gpu_graph.hpp"
#include "cuda_helper.hpp"

struct pool{
    int source; int target; int *d_label; int *d_visited;
    int* d_edges; int* d_dest; int* h_label; int* h_visited; 
    int numVertex; int* h_returned; int *d_frontier;
    int *d_c_frontier_tail; int *d_p_frontier_tail;
    int*check; int*h_check;
    cudaStream_t stream0; cudaStream_t stream1;

    int *d_c_frontier;
    int *d_p_frontier;
    int* state; //for set once
};

template<class Obj>
void BFS_host(pool P,Obj* wrap_obj_graph, gpu_graph_t* _graph); 

__global__ void BFS_Bqueue_kernel(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited); 

__global__ void memset_kernel(int* d_label, int* d_visited, int* d_frontier, int S, int T, int NUM, int* d_p_frontier_tail, int* check); 

__global__ void checking(int num,int* visited, int source, int target, int* check,int *d_c_frontier,int*d_p_frontier);

void run_kernels_graph(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, 
    int* edges, int* dest, int* label, int* visited, 
    int numVertex, int source, int target, int* h_label, int* check, gpu_graph_t &g, cudaStream_t s, int* state);

int degree(int* dest,int* edges,int source); 

int path(int* exclude, int num, int* h_dest, int* h_edges, int state, int elen);

pool init_all(int elen, int dlen, int* devmem);