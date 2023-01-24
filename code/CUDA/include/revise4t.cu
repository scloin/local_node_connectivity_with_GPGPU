#include <revise4t.hpp>
#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 512

void BFS_host(pool P)
{
    int* h_p_frontier_tail=&P.h_visited[2*P.numVertex];
    *h_p_frontier_tail = 2;

    int S = P.source<<1;
    int T = (P.target<<1)+1; 
    int i; int k; int * temp;
    int num_blocks = ((2*P.numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
    memset_kernel<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(P.d_label, P.d_visited, P.d_frontier, S, T, P.numVertex,P.d_p_frontier_tail);
    for (i=0;i<2*P.numVertex;i++) {
        if(P.h_visited[i]!=-1)
            P.h_visited[i] =0;
        }
    P.h_visited[S] =1;P.h_visited[T] =1; 
    CUDA_CHECK(cudaMemcpyAsync(P.d_visited, P.h_visited, 2*P.numVertex*sizeof(int), cudaMemcpyHostToDevice,P.stream)); 

    int *d_c_frontier = &P.d_frontier[0];
    int *d_p_frontier = &P.d_frontier[P.numVertex];
    int check=0; 
    
    while (*h_p_frontier_tail > 0&&check==0) { 

        num_blocks = (*h_p_frontier_tail+BLOCK_SIZE-1) / BLOCK_SIZE;
        BFS_Bqueue_kernel<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(d_p_frontier, P.d_p_frontier_tail,d_c_frontier, P.d_c_frontier_tail, P.d_edges, P.d_dest, P.d_label, P.d_visited);
        //BFS_Bqueue<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(d_p_frontier, P.d_p_frontier_tail,d_c_frontier, P.d_c_frontier_tail, P.d_edges, P.d_dest, P.d_label, P.d_visited);
        
        //CUDA_CHECK(cudaMemcpy(P.h_visited, P.d_visited, (2*P.numVertex+1)*sizeof(int), cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaMemcpyAsync(P.h_label, P.d_label, (4*P.numVertex+1)*sizeof(int), cudaMemcpyDeviceToHost,P.stream)); 

        for(k=0;k<P.numVertex;k++){
            if(P.h_visited[2*k]==1&&P.h_visited[2*k+1]==1&&k!=P.source&&k!=P.target){
                check=1;
                break;
            }
        } 
        temp = d_c_frontier;
        d_c_frontier = d_p_frontier;
        d_p_frontier = temp;


    }
    //CUDA_CHECK(cudaMemcpy(P.h_label, P.d_label, 2*P.numVertex*sizeof(int), cudaMemcpyDeviceToHost)); 

    int min = P.numVertex;
    int meet = -1; 
    //printf("\n...\n");
    for(k=0;k<P.numVertex;k++){// h_label[2*k+1], h_visited[2*k], h_visited[2*k+1]);
        if(P.h_visited[2*k]==1&&P.h_visited[2*k+1]==1&&P.h_label[2*k]*P.h_label[2*k+1]!=0){
            if(min>P.h_label[2*k]+P.h_label[2*k+1]){
                min = P.h_label[2*k]+P.h_label[2*k+1];
                meet =k;
            }
        }
    }
    //if(P.source==0 && P.target==282) printf("meet : %d, visited : [%d,%d], label : [%d, %d]\n", meet, P.h_visited[2*meet], P.h_visited[2*meet+1], P.h_label[2*meet], P.h_label[2*meet+1]);
    *P.h_returned=meet;
} 

__global__ void BFS_Bqueue_kernel(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited){
    __shared__ int c_frontier_s[BLOCK_QUEUE_SIZE];
    __shared__ int c_frontier_tail_s, our_c_frontier_tail; 
    __shared__ int w_queue[16][16], our_c_frontier_s_tail[16];
    __shared__ int w_tail[16];

    if (threadIdx.x <16) {
        c_frontier_tail_s = 0;
        w_tail[threadIdx.x]=0;
    }
    __syncthreads(); 
    const int wid = threadIdx.x%16;
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < *p_frontier_tail) {
        const int my_vertex = p_frontier[tid]>>1;
        const int my_state = p_frontier[tid]&1;
        for (int i = edges[my_vertex]; i < edges[my_vertex+1]; i++) {
            int was_visited=2;
            if(visited[(dest[i]<<1)+my_state]!=-1)
                was_visited = atomicExch(&(visited[(dest[i]<<1)+my_state]), 1);
            if (was_visited==0) {
                label[(dest[i]<<1)+my_state] = label[p_frontier[tid]] + 1;
                //const int my_tail = atomicAdd(&c_frontier_tail_s, 1);
                const int my_w_tail = atomicAdd(&w_tail[wid], 1);
                if (my_w_tail < 16) {
                    w_queue[wid][my_w_tail] = (dest[i]<<1)+my_state;
                }
                else {
                    w_tail[wid]=16;
                    const int my_block_tail = atomicAdd(&c_frontier_tail_s, 1);
                    if (my_block_tail < BLOCK_QUEUE_SIZE) {
                        c_frontier_s[my_block_tail] = (dest[i]<<1)+my_state;
                    }
                    else {
                        c_frontier_tail_s = BLOCK_QUEUE_SIZE;
                        const int my_global_tail = atomicAdd(c_frontier_tail, 1);
                        c_frontier[my_global_tail] = (dest[i]<<1)+my_state;
                    }
                    //c_frontier_s[my_block_tail] = (dest[i]<<1)+my_state;
                }
            }
        }
    }
    __syncthreads(); 

    if (threadIdx.x <16) {
        our_c_frontier_s_tail[wid] = atomicAdd(&c_frontier_tail_s, w_tail[wid]);
        //printf("our_c_- %d, w_tail%d\n", our_c_frontier_s_tail[wid], w_tail[wid]);
    }
    __syncthreads(); 

    for (int i = threadIdx.x/16; i < w_tail[wid]; i += 16) {
        //printf("tid=%d, wqueue[%d][%d]=%d, c_frontier_s[%d]\n", threadIdx.x,wid,i,w_queue[wid][i],our_c_frontier_s_tail[wid] + i);
        if(our_c_frontier_s_tail[wid] + i>=BLOCK_QUEUE_SIZE) {
            const int second_global_tail = atomicAdd(c_frontier_tail, 1);
            c_frontier[second_global_tail] = w_queue[wid][i];
            }
        else{
            c_frontier_s[our_c_frontier_s_tail[wid] + i] = w_queue[wid][i];
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
    __syncthreads();
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

pool init_pool(int elen, int dlen, int* devmem, int* d_edges){

    pool P;
    CUDA_CHECK( cudaStreamCreate(&P.stream) );
    P.numVertex=elen-1;

    P.d_edges=d_edges;
    P.d_dest=&d_edges[elen];

    P.h_label=(int*)malloc((4*(P.numVertex)+2)*sizeof(int)); 
    //CUDA_CHECK(cudaMallocHost((void**)&P.h_label, (4*(P.numVertex)+2)*sizeof(int)));
    P.h_visited=&(P.h_label[2*P.numVertex]);

    P.d_label=devmem;
    P.d_visited=&(devmem[2*P.numVertex]);

    P.d_frontier=&(devmem[4*P.numVertex+1]);
    P.d_c_frontier_tail=&(devmem[6*P.numVertex+1]);
    P.d_p_frontier_tail=&(devmem[4*P.numVertex]);

    P.h_returned=(int*)malloc(sizeof(int));
    P.source=0;
    P.target=1; 

    return P;
}

void compute(int* h_dest,int * h_edges, pool P,FILE* fp1){
    int *exclude_S=(int *)malloc(P.numVertex*sizeof(int));
    int *exclude_T=(int *)malloc(P.numVertex*sizeof(int)); 
    int init = P.target; 
<<<<<<< HEAD
    int count; int N; int* Ad; int num; int i; int j; int KK;
    while(P.source<P.numVertex){
        P.target=P.source+init;
        while(P.target<P.numVertex){
    count=0; 
=======
    while(P.source<P.numVertex){
        P.target=P.source+init;
        while(P.target<P.numVertex){
    int count=0; 
>>>>>>> beeae6054c79a5b3d5057253c7de5bd990b5e31f
    //verify that it is connected directly
    N = degree(h_dest,h_edges,P.source);
    Ad=(int*)malloc(N*sizeof(int));
    memcpy(Ad,&(h_dest[h_edges[P.source]]),N*sizeof(int)); 

    //int check=0;
    for(j = 0; j<N;j++){
        if(Ad[j]==P.target){
            P.h_visited[(P.source<<1)+1]=-1;
            P.h_visited[(P.target<<1)]=-1;
            count++; 
            break;
        }
    }
    free(Ad);
    //loop while count every component
    while(1){ 
    
<<<<<<< HEAD
        //BFS_host(P);
        BFS_host2(P);
        num = *P.h_returned;
=======
        BFS_host(P);
        //std::list<int> list;//
        int num = *P.h_returned;
>>>>>>> beeae6054c79a5b3d5057253c7de5bd990b5e31f
        if(num==-1) {
            break;
            }
        count++;
<<<<<<< HEAD

        i=num;
        j=num; 
=======
        //std::list<int>::iterator begin_iter = list.begin();//
        //std::list<int>::iterator end_iter = list.end();//
        //list.insert(end_iter, num);//
        //begin_iter--; //

        int i=num; 
>>>>>>> beeae6054c79a5b3d5057253c7de5bd990b5e31f
        P.h_visited[i<<1]=-1;
        P.h_visited[(i<<1)+1]=-1;

        for (KK=0; KK<P.numVertex;KK++){
                exclude_S[KK] = P.h_label[KK*2];
                exclude_T[KK] = P.h_label[KK*2+1];
            } 

<<<<<<< HEAD
        int tempi=0;
        int tempj=0;
        int error =0;
        while((i>-1)||(j>-1)){
            if(i>-1){
                if(i!=num) {
                P.h_visited[i<<1]=-1;
                P.h_visited[(i<<1)+1]=-1;
                }
                tempi=i;
                i = path(exclude_T,i,h_dest,h_edges,0,P.numVertex);
                if(tempi==i){
                    error=1;
                    break;
                }
            }
            if(j>-1){
                if(j!=num) {
                P.h_visited[j<<1]=-1;
                P.h_visited[(j<<1)+1]=-1;
                }
                tempj=j;
                j = path(exclude_S,j,h_dest,h_edges,1,P.numVertex);
                if(tempj==j){
                    error=1;
                    break;
                }
            }
        }
        if(error==1){
            count--;
            break;
        }
    }
    //record count
    fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    
    P.target+=4;
    //P.target++;
    for (i=0;i<2*P.numVertex;i++) {
        P.h_visited[i] =0;
        }
        }
    P.source++;
    }

    free(exclude_S);
    free(exclude_T);
}

__global__ void BFS_Bqueue(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited){
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
=======
        while((i>-1)){
            if(i!=num) {
            P.h_visited[i<<1]=-1;
            P.h_visited[(i<<1)+1]=-1;
            }
            i = path(exclude_T,i,h_dest,h_edges,0,P.numVertex);
            } 
        i=num;
        while((i>-1)){
            if(i!=num) {
            P.h_visited[i<<1]=-1;
            P.h_visited[(i<<1)+1]=-1;
            }
            i = path(exclude_S,i,h_dest,h_edges,1,P.numVertex); 
>>>>>>> beeae6054c79a5b3d5057253c7de5bd990b5e31f

/*
BFS with Cuda c++, using shared memory, and warp level queue with less atomic operation
*/
__global__ void BFS_less_atomic(int* d_edges, int* d_dest, int* d_label, int* d_visited, int numVertex, int level){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    const int n = 2*numVertex;
    if (tid < n) {
        const int my_vertex = tid>>1;
        const int my_state = tid&1;
        if (d_label[tid] == level) {
            for (int i = d_edges[my_vertex]; i < d_edges[my_vertex+1]; i++) {
                int v = (d_dest[i]<<1)+my_state;
                if (d_visited[v] == 0) {
                    d_label[v] = level + 1;
                    d_visited[v] = 1;
                }
            }
        }
    }
}

/*
Revise of BFS_host:
    - delete the frontier queue
    - use BFS_less_atomic() to replace BFS_Bqueue_kernel()
*/
void BFS_host2(pool P)
{
    int S = P.source<<1;
    int T = (P.target<<1)+1; 
    int i; int k;
    int num_blocks = ((2*P.numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
    memset_kernel2<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(P.d_label, P.d_visited, S, T, P.numVertex);

    //init and copy visited
    for (i=0;i<2*P.numVertex;i++) {
        if(P.h_visited[i]!=-1)
            P.h_visited[i] =0;
        }
    P.h_visited[S] =1;P.h_visited[T] =1; 
    CUDA_CHECK(cudaMemcpyAsync(P.d_visited, P.h_visited, 2*P.numVertex*sizeof(int), cudaMemcpyHostToDevice,P.stream)); 

    int check=0; 
    int level = 0;

    num_blocks = ((2*P.numVertex)+BLOCK_SIZE-1) / BLOCK_SIZE;
    int csum =0;
    while (check==0) { 
        //cudaMemcpyAsync(d_done, &true_value, sizeof(int), cudaMemcpyHostToDevice,P.stream);
        //CUDA_CHECK(cudaMemsetAsync(d_done, 0, sizeof(int),P.stream));

        BFS_less_atomic<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(P.d_edges, P.d_dest, P.d_label, P.d_visited, P.numVertex, level);
        cudaDeviceSynchronize();
        CUDA_CHECK(cudaMemcpyAsync(P.h_label, P.d_label, (4*P.numVertex+1)*sizeof(int), cudaMemcpyDeviceToHost,P.stream)); 
        
        cudaDeviceSynchronize();
        for(k=0;k<P.numVertex;k++){
            if ((P.h_label[2*k]==level)||(P.h_label[2*k+1]==level)){
                csum++;
            };
            if(P.h_visited[2*k]==1&&P.h_visited[2*k+1]==1&&k!=P.source&&k!=P.target){
                check=1;
                break;
            }
<<<<<<< HEAD
        } 
        if(csum==0){
            check=1;
        }
        csum=0;
        level++;
    }

    int min = P.numVertex;
    int meet = -1; 
    for(k=0;k<P.numVertex;k++){
        if(P.h_visited[2*k]==1&&P.h_visited[2*k+1]==1&&P.h_label[2*k]*P.h_label[2*k+1]!=0){
            if(min>P.h_label[2*k]+P.h_label[2*k+1]){
                min = P.h_label[2*k]+P.h_label[2*k+1];
                meet =k;
            }
        }
    }

    *P.h_returned=meet;
} 

__global__ void memset_kernel2(int* d_label, int* d_visited, int S, int T, int NUM){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < NUM*2) {
        *(d_label+tid) = -1;
        if (tid==S||tid==T) {
            *(d_label+tid) = 0;
                }
    }

=======
    }
    //record count
    fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    
    P.target+=3;
    for (int i=0;i<2*P.numVertex;i++) {
        P.h_visited[i] =0;
        }
        }
    P.source++;
    }

    free(exclude_S);
    free(exclude_T);
>>>>>>> beeae6054c79a5b3d5057253c7de5bd990b5e31f
}