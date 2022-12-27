/*
add POOL
merge memcpy
[X] bfsmain으로 loop end check를 kernel로
loopend 를 kernel로
graph 적용 > done!
overwrap test..
*/

#include "revise5a.hpp"

#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 32


using namespace std;

int main(){

    FILE* fp = fopen("./graph/graph.txt","r");
    int dlen;
    int elen; 
    int * h_dest;
    int * h_edges; 

    fscanf(fp, "%d ", &dlen);
    fscanf(fp, "%d ", &elen); 

    gpu_graph_t _graph;
    gpu_graph_t _graph1;

    int *h_data = (int*)malloc((elen+dlen)*sizeof(int)); 
    h_edges= h_data;
    h_dest= &h_data[elen];
    int k =0;
    int i;
    while(k<dlen){
        fscanf(fp, "%d ", &i);
        h_dest[k]=i;
        k++;
    }
    k=0;
    while(k<elen){
        fscanf(fp, "%d ", &i);
        h_edges[k]=i;
        k++;
    }
    fclose(fp);

    int* devmem;
    CUDA_CHECK(cudaMalloc((void**)&devmem, (elen+dlen+6*(elen-1)+3)*sizeof(int)));

    int* devmem1;
    CUDA_CHECK(cudaMalloc((void**)&devmem1, (elen+dlen+6*(elen-1)+3)*sizeof(int)));

    pool P=init_all(elen, dlen, devmem);
    pool P1=init_all(elen, dlen, devmem1);

    CUDA_CHECK(cudaMemcpy(P1.d_edges, h_edges, (dlen+elen)*sizeof(int), cudaMemcpyHostToDevice)); 

    P.d_edges = P1.d_edges;
    P.d_dest = P1.d_dest;

    int *exclude_S=(int *)malloc(P.numVertex*sizeof(int));
    int *exclude_T=(int *)malloc(P.numVertex*sizeof(int)); 

    FILE* fp1 = fopen("result/gtest.txt","w"); 

    auto wrap_obj_graph = [&](gpu_graph_t &g, cudaStream_t s) {
        pool C = P;
        run_kernels_graph(C.d_p_frontier, C.d_p_frontier_tail, C.d_c_frontier, C.d_c_frontier_tail,
            C.d_edges, C.d_dest, C.d_label,C.d_visited,
            C.numVertex, C.source, C.target, C.h_label, C.check, g, s,C.state);
    };

    auto wrap_obj_graph1 = [&](gpu_graph_t &g, cudaStream_t s) {
        pool C = P1;
        run_kernels_graph(C.d_p_frontier, C.d_p_frontier_tail, C.d_c_frontier, C.d_c_frontier_tail,
            C.d_edges, C.d_dest, C.d_label,C.d_visited,
            C.numVertex, C.source, C.target, C.h_label, C.check, g, s,C.state);
    };

    int* Ad;

    while(P.source<P.numVertex){
        P.target=P.source+1;
        P1.source=P.source;
        P1.target=P.target;

    while(P.target<P.numVertex){
    for (int i=0;i<2*P.numVertex;i++) {
        P.h_visited[i] =0;
    }
    int count=0; 

    int N = degree(h_dest,h_edges,P.source);

    Ad=&(h_dest[h_edges[P.source]]);

    int check=0;
    for(int j = 0; j<N;j++){
        if(*(Ad+j)==P.target){
            check=1; 

            break;
        }
    }
    if(check){
        P.h_visited[(P.source<<1)+1]=-1;
        P.h_visited[(P.target<<1)]=-1;
        count++;
    }
    //thread t1;
    thread t2;
    while(1){ 
    BFS_host(P, &wrap_obj_graph, &_graph);

    //thread test.. 
    t2=thread{degree,h_dest,h_edges,P.source};
    //


    //t1.join(); 
    t2.join();
    list<int> list1;
    int num = *P.h_returned;
    if(num==-1) {
        break;
        }
    count++;
    list<int>::iterator begin_iter = list1.begin();
    list<int>::iterator end_iter = list1.end();
    list1.insert(end_iter, num);
    begin_iter--; 

    int i=num; 

    for (int KK=0; KK<P.numVertex;KK++){
            exclude_S[KK] = P.h_label[KK*2];
            exclude_T[KK] = P.h_label[KK*2+1];
        } 

    while((i>-1)){
        if(i!=num) {
        list1.insert(end_iter, i);}

        i = path(exclude_T,i,h_dest,h_edges,0,elen-1);
        } 

    i=num;
    while((i>-1)){
        if(i!=num) {
        list1.insert(begin_iter, i);
        begin_iter--;}
        i = path(exclude_S,i,h_dest,h_edges,1,elen-1);
        }
    while (list1.empty()==0) {
        i=list1.front();
        P.h_visited[i<<1]=-1;
        P.h_visited[(i<<1)+1]=-1;
        list1.pop_front(); 

    } 

    list1.clear();

    }
    
    
    fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    P.target++;

    } 

    P.source++;}
    fclose(fp1); 
    
    _graph.dest_graph();
    //free

    CUDA_CHECK(cudaFree(devmem));
    CUDA_CHECK(cudaFree(devmem1));
    free(h_data);

    free(P.h_label);
    free(P1.h_label);

    free(P.h_returned);
    free(P1.h_returned);
    free(exclude_T);
    free(exclude_S);
    CUDA_CHECK( cudaStreamDestroy(P.stream0));
    CUDA_CHECK( cudaStreamDestroy(P1.stream0));
    cudaDeviceReset();
    return 0;
} 


template<class Obj>
void BFS_host(pool P, Obj* wrap_obj_graph,gpu_graph_t* _graph)
{

    int* h_p_frontier_tail=&P.h_label[4*P.numVertex];
    *h_p_frontier_tail = 2;

    int S = P.source<<1;
    int T = (P.target<<1)+1; 

    int num_blocks = ((2*P.numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
    
    for (int i=0;i<2*P.numVertex;i++) {
        if(P.h_visited[i]!=-1)
            P.h_visited[i]=0;
        }
    P.h_visited[S] =1;P.h_visited[T] =1; 
    CUDA_CHECK(cudaMemcpy(P.d_visited, P.h_visited, 2*P.numVertex*sizeof(int), cudaMemcpyHostToDevice)); 
    memset_kernel<<<num_blocks, BLOCK_SIZE,0,P.stream0>>>(P.d_label, P.d_visited, P.d_frontier, S, T, P.numVertex,P.d_p_frontier_tail, P.check);

    *P.h_check=0; 

    while (*h_p_frontier_tail > 0&&*P.h_check==0) { 
        (*_graph).wrap(*wrap_obj_graph, P.stream0);
        cudaDeviceSynchronize();
    }
    int min = P.numVertex;
    int meet = -1; 
    *P.state=2;
    for(int k=0;k<P.numVertex;k++){
        if(P.h_visited[2*k]==1&&P.h_visited[2*k+1]==1&&k!=P.source&&k!=P.target){
            if(min>P.h_label[2*k]+P.h_label[2*k+1]&&(P.h_label[2*k]*P.h_label[2*k+1])>0){
                min = P.h_label[2*k]+P.h_label[2*k+1];
                meet =k;
            }
        }
    }
    //printf("====>>%d\n", min);
    *P.h_returned=meet;
} 

__global__ void BFS_Bqueue_kernel(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited){
    __shared__ int c_frontier_s[BLOCK_QUEUE_SIZE];
    __shared__ int c_frontier_tail_s, our_c_frontier_tail; 

    if (threadIdx.x == 0)
        c_frontier_tail_s = 0;
    __syncthreads(); 

    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < *p_frontier_tail) {
        const int my_vertex = p_frontier[tid]>>1;
        const int my_state = p_frontier[tid]&1;
        for (int i = edges[my_vertex]; i < edges[my_vertex+1]; i++) {
            int was_visited=2;
            if(visited[(dest[i]<<1)+my_state]!=-1)
                was_visited = atomicExch(&(visited[(dest[i]<<1)+my_state]), 1);
            if (!was_visited) {
                label[(dest[i]<<1)+my_state] = label[p_frontier[tid]] + 1;
                const int my_tail = atomicAdd(&c_frontier_tail_s, 1);
                if (my_tail < BLOCK_QUEUE_SIZE) {
                    c_frontier_s[my_tail] = (dest[i]<<1)+my_state;
                }
                else {
                    c_frontier_tail_s = BLOCK_QUEUE_SIZE;
                    const int my_global_tail = atomicAdd(c_frontier_tail, 1);
                    c_frontier[my_global_tail] = (dest[i]<<1)+my_state;
                }
            }
        }
    }
    __syncthreads(); 

    if (threadIdx.x == 0) {
        our_c_frontier_tail = atomicAdd(c_frontier_tail, c_frontier_tail_s);
    }
    __syncthreads(); 

    for (int i = threadIdx.x; i < c_frontier_tail_s; i += blockDim.x) {
        c_frontier[our_c_frontier_tail + i] = c_frontier_s[i];
    }
    //__syncthreads();
        if (tid == 0) {
        *p_frontier_tail = atomicExch(c_frontier_tail, 0);
    }
} 

int path(int* exclude, int num, int* h_dest, int* h_edges, int state,int elen){
    int numT=num;
    int check=0;
    int i=0; 

    int N = degree(h_dest,h_edges,num); 

    if(exclude[numT]==1){
        if(state==0) {
            return -2;
        }else {
            return -1;
        }
    }
    else{
    int* AA=(int*)malloc(N*sizeof(int));
    memcpy(AA,&(h_dest[h_edges[num]]),N*sizeof(int));
        for (i = 0; i < elen; i++){
            if (exclude[i] == exclude[numT]-1) {
                for(int j = 0; j<degree(h_dest,h_edges,numT);j++){
                    if(AA[j]==i){
                        check=1;
                        break;
                    }
                }
            }
            if(check==1) break;
        }
    free(AA);
    }
    if(check!=1) i=numT; 

    return i;
} 

int degree(int* dest,int* edges,int source){ 

    return edges[source+1]-edges[source];
} 

__global__ void memset_kernel(int* d_label, int* d_visited, int* d_frontier, int S, int T, int NUM, int* d_p_frontier_tail, int * check){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < NUM*2) {
        *(d_label+tid) = -1;
        //*(d_visited+tid) = 0;
        *(d_frontier+tid) =0;
        if (tid==S||tid==T) {
            *(d_label+tid) = 0;
                }
        if (tid ==NUM) *(d_frontier+tid) = S;
        if (tid ==NUM+1) *(d_frontier+tid) = T;
                *(d_p_frontier_tail) =2;
    }
    *check=0;

}

__global__ void checking(int num,int* visited, int source, int target, int* check,int *d_c_frontier,int*d_p_frontier){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if(tid<num)
        if(visited[2*tid]==1&&visited[2*tid+1]==1&&tid!=source&&tid!=target
){

            *check=1;
        }

}

/* 
1. 파라미터 변경필요 V
2. memcpy node 추가 V
3. statci kernel 추가 V
>>static도 graph에 포함되는지 확인필요!
>>>실행해보기!!!
*/
void run_kernels_graph(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail,  
    int* edges, int* dest, int* label, int* visited, 
    int numVertex, int source, int target, int* h_label, int* check, gpu_graph_t &g, cudaStream_t s, int* state){
    //constexpr int threads = 256;
    //int blocks = (size + threads - 1) / threads;
    int num_blocks = (h_label[4*numVertex]+BLOCK_SIZE-1) / BLOCK_SIZE;
    int num_blocks1 = ((numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
        
    cudaKernelNodeParams params;
        params.blockDim = {static_cast<unsigned int>(BLOCK_SIZE), 1, 1};
        params.gridDim = {static_cast<unsigned int>(num_blocks), 1, 1};
        params.sharedMemBytes = 0;
        params.func = reinterpret_cast<void *>(BFS_Bqueue_kernel);
        void *args[] = {&p_frontier, &p_frontier_tail, &c_frontier, &c_frontier_tail, &edges, &dest, &label, &visited};
        params.kernelParams = args;
        params.extra = nullptr;

    cudaKernelNodeParams params1;
        params1.blockDim = {static_cast<unsigned int>(BLOCK_SIZE), 1, 1};
        params1.gridDim = {static_cast<unsigned int>(num_blocks1), 1, 1};
        params1.sharedMemBytes = 0;
        params1.func = reinterpret_cast<void *>(checking);
        void *args1[] = {&numVertex, &visited, &source, &target, &check, &c_frontier ,&p_frontier};
        params1.kernelParams = args1;
        params1.extra = nullptr;

    void* temp;
    //if (g.state() == gpu_graph_t::state_t::capture) {
    if (*state==0) {
      g.add_kernel_node(0, params, s);
      g.add_kernel_node(1, params1, s);
      g.add_memcpy_node(2, label, h_label, (4*numVertex+2)*sizeof(int), s);
      *state=1;
    //} else if (g.state() == gpu_graph_t::state_t::update) {
    } else if (*state==1) {
        temp = params.kernelParams[0];
        params.kernelParams[0]= params.kernelParams[2];
        params.kernelParams[2]=temp;

      g.update_kernel_node(0, params);
      g.update_kernel_node(1, params1);
      *state=2;
    }
    else if (*state==2) {
      g.update_kernel_node(0, params);
      g.update_kernel_node(1, params1);
      *state=1;
    }

}

pool init_all(int elen, int dlen, int* devmem){
    pool P;
    CUDA_CHECK( cudaStreamCreate(&P.stream0) );
    P.numVertex=elen-1;
    P.h_label=(int*)malloc((4*(P.numVertex)+2)*sizeof(int)); 
    P.h_visited=&(P.h_label[2*P.numVertex]);
    P.h_check  =&(P.h_label[4*P.numVertex+1]);
    P.state=(int*)malloc(sizeof(int)); 
    *P.state=0;

    P.d_edges=devmem;
    P.d_dest=&(devmem[elen]);
    P.d_label=&(devmem[dlen+elen]);
    P.d_visited=&(P.d_label[2*P.numVertex]);

    P.d_p_frontier_tail=&(P.d_label[4*P.numVertex]);
    P.check=&(P.d_label[4*P.numVertex+1]);

    P.d_frontier=&(P.d_label[4*P.numVertex+2]);
    P.d_c_frontier_tail=&(P.d_label[6*P.numVertex+2]);

    P.d_c_frontier = &P.d_frontier[0];
    P.d_p_frontier = &P.d_frontier[P.numVertex];

    P.h_returned=(int*)malloc(sizeof(int));
    P.source=0;
    P.target=1; 

    return P;
}