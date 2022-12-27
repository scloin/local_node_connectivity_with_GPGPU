#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <list>
#include <common/common.h>
#include <common/common_string.h> 

#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 32 

void BFS_host(int source, int target, int *d_label, int *d_visited, int* d_edges, int* d_dest, int* h_label,int* h_visited, int numVertex, int* h_returned, int *d_frontier,int *d_c_frontier_tail,int *d_p_frontier_tail ,cudaStream_t stream0 ,cudaStream_t stream1); 

__global__ void BFS_Bqueue_kernel(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited); 

__global__ void memset_kernel(int* d_label, int* d_visited, int* d_frontier, int S, int T, int NUM, int* d_p_frontier_tail); 

int degree(int* dest,int* edges,int source); 

int path(int* exclude, int num, int* h_dest, int* h_edges, int state, int elen); 

int main(){
    cudaStream_t stream0; CUDA_CHECK( cudaStreamCreate(&stream0) );
    cudaStream_t stream1; CUDA_CHECK( cudaStreamCreate(&stream1) );
    FILE* fp = fopen("./graph/graph.txt","r");
    int dlen;
    int elen; 

    int * h_dest;
    int * h_edges; 

    fscanf(fp, "%d ", &dlen);
    fscanf(fp, "%d ", &elen); 

    h_dest = (int*)malloc(dlen*sizeof(int));
    h_edges = (int*)malloc(elen*sizeof(int)); 

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
        //printf("%d\n",i);
        h_edges[k]=i;
        k++;
    }
    fclose(fp);
    int *h_label;
    h_label=(int*)malloc(2*(elen-1)*sizeof(int)); 

    int *h_visited;
    h_visited=(int*)malloc(2*(elen-1)*sizeof(int)); 

    int *d_edges, *d_dest, *d_label, *d_visited;
    CUDA_CHECK(cudaMalloc((void**)&d_edges, (elen)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_dest, dlen*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_label, 2*(elen-1)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_visited, 2*(elen-1)*sizeof(int)));


    int *d_frontier, *d_c_frontier_tail, *d_p_frontier_tail;
    CUDA_CHECK(cudaMalloc((void**)&d_frontier, 2*(elen-1)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_c_frontier_tail, sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&d_p_frontier_tail, sizeof(int))); 

    CUDA_CHECK(cudaMemcpyAsync(d_edges, h_edges, (elen)*sizeof(int), cudaMemcpyHostToDevice,stream0)); 

    CUDA_CHECK(cudaMemcpyAsync(d_dest, h_dest, dlen*sizeof(int), cudaMemcpyHostToDevice,stream1)); 

    int *exclude_S=(int *)malloc((elen-1)*sizeof(int));
    int *exclude_T=(int *)malloc((elen-1)*sizeof(int)); 

    int * h_returned=(int*)malloc(sizeof(int));
    int source=0;
    int target=0; 

    FILE* fp1 = fopen("result/test2.txt","w"); 

    while(source<(elen-1)){
    target=source+1;
    while(target<(elen-1)){
    int count=0; 
    int N = degree(h_dest,h_edges,source);
    int* Ad=(int*)malloc(N*sizeof(int));
    memcpy(Ad,&(h_dest[h_edges[source]]),N*sizeof(int)); 

    int check=0;
    for(int j = 0; j<N;j++){
        if(Ad[j]==target){
            check=1; 

            break;
        }
    }
    if(check==1){
        h_visited[(source<<1)+1]=-1;
        h_visited[(target<<1)]=-1;
        count++;
    }
    free(Ad);
    if(1) {
    

    while(1){ 

    BFS_host(source, target, d_label, d_visited, d_edges, d_dest, h_label, h_visited, (elen-1), h_returned, d_frontier, d_c_frontier_tail, d_p_frontier_tail,stream0,stream1); 

    std::list<int> list1;
    int num = *h_returned;
    if(num==-1) {
        break;
        }
    count++;
    std::list<int>::iterator begin_iter = list1.begin();
    std::list<int>::iterator end_iter = list1.end();
    list1.insert(end_iter, num);
    begin_iter--; 

    int i=num; 

    for (int KK=0; KK<(elen-1);KK++){
            exclude_S[KK] = h_label[KK*2];
            exclude_T[KK] = h_label[KK*2+1];
        } 

    while((i>-1)){
        if(i!=num) {
        list1.insert(end_iter, i);}
        i = path(exclude_T,i,h_dest,h_edges,0,elen-1);
        } 

    i=num;
        for (int i=0; i<(elen-1);i++)
            printf("%2d ", exclude_S[i]);
    printf("\n");
    while((i>-1)){
        if(i!=num) {
        list1.insert(begin_iter, i);
        begin_iter--;}
        i = path(exclude_S,i,h_dest,h_edges,1,elen-1); 

        }
    while (list1.empty()==0) {
        i=list1.front();
        h_visited[i<<1]=-1;
        h_visited[(i<<1)+1]=-1;
        list1.pop_front(); 

    } 

    list1.clear();
    }//무한 while
    //check문
    }

    fprintf(fp1,"[%d, %d] %d\n", source, target, count);
    target++;

    // for (int i=0; i<2*(elen-1);i++)
    //     printf("%2d ", h_label[i]);
    // printf("\n");

    for (int i=0;i<2*(elen-1);i++) {
        h_visited[i] =0;
        }
    } 

    source++;}
    fclose(fp1); 

    //free
    CUDA_CHECK(cudaFree(d_edges));
    CUDA_CHECK(cudaFree(d_dest));
    CUDA_CHECK(cudaFree(d_label));
    CUDA_CHECK(cudaFree(d_visited)); 

    CUDA_CHECK(cudaFree(d_frontier));
    CUDA_CHECK(cudaFree(d_c_frontier_tail));
    CUDA_CHECK(cudaFree(d_p_frontier_tail)); 

    free(h_dest);
    free(h_edges);
    cudaFreeHost(h_label);
    cudaFreeHost(h_visited);
    free(h_returned);
    free(exclude_T);
    free(exclude_S);
    CUDA_CHECK( cudaStreamDestroy(stream0));
    CUDA_CHECK( cudaStreamDestroy(stream1));
    cudaDeviceReset();
    return 0;
} 

void BFS_host(int source, int target, int *d_label, int *d_visited, int* d_edges, int* d_dest, int* h_label,int* h_visited, int numVertex, int* h_returned, int *d_frontier,int *d_c_frontier_tail,int *d_p_frontier_tail,cudaStream_t stream0 ,cudaStream_t stream1)
{
    int h_p_frontier_tail = 2; 

    int S = source<<1;
    int T = (target<<1)+1; 

    int num_blocks = ((2*numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
    memset_kernel<<<num_blocks, BLOCK_SIZE, 0, stream0>>>(d_label, d_visited, d_frontier, S, T, numVertex,d_p_frontier_tail);
    for (int i=0;i<2*numVertex;i++) {
        if(h_visited[i]!=-1)
            h_visited[i] =0;
        }
    h_visited[S] =1;h_visited[T] =1; 
    CUDA_CHECK(cudaMemcpyAsync(d_visited, h_visited, 2*numVertex*sizeof(int), cudaMemcpyHostToDevice,stream1)); 

    int *d_c_frontier = &d_frontier[0];
    int *d_p_frontier = &d_frontier[numVertex];
    int check=0; 

    while (h_p_frontier_tail > 0&&check==0) { 

        num_blocks = (h_p_frontier_tail+BLOCK_SIZE-1) / BLOCK_SIZE;
        BFS_Bqueue_kernel<<<num_blocks, BLOCK_SIZE>>>(d_p_frontier, d_p_frontier_tail, d_c_frontier, d_c_frontier_tail, d_edges, d_dest, d_label, d_visited);
        CUDA_CHECK(cudaMemcpyAsync(&h_p_frontier_tail, d_c_frontier_tail, sizeof(int), cudaMemcpyDeviceToHost,stream0));
        CUDA_CHECK(cudaMemcpyAsync(d_p_frontier_tail, d_c_frontier_tail, sizeof(int), cudaMemcpyDeviceToDevice,stream1));
        CUDA_CHECK(cudaMemcpyAsync(h_visited, d_visited, 2*numVertex*sizeof(int), cudaMemcpyDeviceToHost,stream0));
        CUDA_CHECK(cudaMemsetAsync(d_c_frontier_tail, 0, sizeof(int),stream1)); 

        for(int k=0;k<numVertex;k++){
            if(h_visited[2*k]==1&&h_visited[2*k+1]==1&&k!=source&&k!=target){
                check=1;
                break;
            }
        } 

        int* temp = d_c_frontier;
        d_c_frontier = d_p_frontier;
        d_p_frontier = temp;


    }
    CUDA_CHECK(cudaMemcpyAsync(h_label, d_label, 2*numVertex*sizeof(int), cudaMemcpyDeviceToHost, stream0)); 

    int min = numVertex;
    int meet = -1; 

    for(int k=0;k<numVertex;k++){// h_label[2*k+1], h_visited[2*k], h_visited[2*k+1]);
        if(h_visited[2*k]==1&&h_visited[2*k+1]==1){
            if(min>h_label[2*k]+h_label[2*k+1]){
                min = h_label[2*k]+h_label[2*k+1];
                meet =k;
            }
        }
    }
    *h_returned=meet;
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

__global__ void memset_kernel(int* d_label, int* d_visited, int* d_frontier, int S, int T, int NUM, int* d_p_frontier_tail){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < NUM*2) {
        *(d_label+tid) = -1;
        *(d_frontier+tid) =0;
        if (tid==S||tid==T) {
            *(d_label+tid) = 0;
                }
        if (tid ==NUM) *(d_frontier+tid) = S;
        if (tid ==NUM+1) *(d_frontier+tid) = T;
                *(d_p_frontier_tail) =2;
    }

}
