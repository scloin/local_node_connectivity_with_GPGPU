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
    memset_kernel<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(S, T, P.d_label, P.d_frontier, P.numVertex, P.d_p_frontier_tail);

    for (i=0;i<2*P.numVertex;i++) {
        if(P.h_visited[i]!=-1)
            P.h_visited[i] =0;
        }

    P.h_visited[P.source] =1;P.h_visited1[P.target] =1; 

    CUDA_CHECK(cudaMemcpyAsync(P.d_visited, P.h_visited, 2*P.numVertex*sizeof(int), cudaMemcpyHostToDevice,P.stream)); 

    int *d_c_frontier = &P.d_frontier[0];
    int *d_p_frontier = &P.d_frontier[P.numVertex];
    int check=0; 
    while (*h_p_frontier_tail > 0&&check==0) { 

        num_blocks = (*h_p_frontier_tail+BLOCK_SIZE-1) / BLOCK_SIZE;

        BFS_Bqueue<<<num_blocks, BLOCK_SIZE,0,P.stream>>>(d_p_frontier, P.d_p_frontier_tail,d_c_frontier, P.d_c_frontier_tail, P.d_edges, P.d_dest, P.d_label, P.d_visited, P.d_label1, P.d_visited1);

        CUDA_CHECK(cudaMemcpyAsync(P.h_label, P.d_label, (4*P.numVertex+1)*sizeof(int), cudaMemcpyDeviceToHost,P.stream)); 
        cudaDeviceSynchronize();



        for(k=0;k<P.numVertex;k++){
            if(P.h_visited[k]==1&&P.h_visited1[k]==1&&k!=P.source&&k!=P.target){
                check=1;
                // printf("\n");
                // for(k=0;k<P.numVertex;k++){
                //     printf("%2d ", P.h_label[k]);
                // }
                //     printf("\n");
                // for(k=0;k<P.numVertex;k++){
                //     printf("%2d ", P.h_label1[k]);
                // }
                // printf("\n");

                // for(k=0;k<P.numVertex;k++){
                //     printf("%2d ", P.h_visited[k]);
                // }
                //     printf("\n");
                // for(k=0;k<P.numVertex;k++){
                //     printf("%2d ", P.h_visited1[k]);
                // }
                // printf("\n");

                break;
            }
        } 
        temp = d_c_frontier;
        d_c_frontier = d_p_frontier;
        d_p_frontier = temp;


    }


    int min = P.numVertex;
    int meet = -1; 
    //printf("\n...\n");
    for(k=0;k<P.numVertex;k++){// h_label[2*k+1], h_visited[2*k], h_visited[2*k+1]);
        if(P.h_visited[k]==1&&P.h_visited1[k]==1&&P.h_label[k]*P.h_label1[k]!=0){
            if(min>P.h_label[k]+P.h_label1[k]){
                min = P.h_label[k]+P.h_label1[k];
                meet =k;
            }
        }
    }

    *P.h_returned=meet;
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

__global__ void memset_kernel(int S, int T, int* d_label,int* d_frontier,int numVertex,int* d_p_frontier_tail){

    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < numVertex*2) {
        *(d_label+tid) = -1;
        *(d_frontier+tid) =0;
        if (tid==(S>>1)||tid==(T>>1)+numVertex) {
            *(d_label+tid) = 0;
                }
        if (tid ==numVertex) *(d_frontier+tid) = S;
        if (tid ==numVertex+1) *(d_frontier+tid) = T;
                *(d_p_frontier_tail) =2;
    }

}

pool init_pool(int elen, int dlen, int* devmem, int* d_edges){

    pool P;
    CUDA_CHECK( cudaStreamCreate(&P.stream) );
    P.numVertex=elen-1;

    P.d_edges=d_edges;
    P.d_dest=&d_edges[elen];


    CUDA_CHECK(cudaMallocHost((void**)&P.h_label, (4*(P.numVertex)+2)*sizeof(int)));
    P.h_label1 =&(P.h_label[P.numVertex]);
    P.h_visited=&(P.h_label[2*P.numVertex]);
    P.h_visited1=&(P.h_label[3*P.numVertex]);

    P.d_label=devmem;
    P.d_label1=&(devmem[P.numVertex]);
    P.d_visited=&(devmem[2*P.numVertex]);
    P.d_visited1=&(devmem[3*P.numVertex]);

    P.d_frontier=&(devmem[4*P.numVertex+1]);
    P.d_c_frontier_tail=&(devmem[6*P.numVertex+1]);
    P.d_p_frontier_tail=&(devmem[4*P.numVertex]);

    P.h_returned=(int*)malloc(sizeof(int));
    P.source=0;
    P.target=1; 

    return P;
}

void compute(int* h_dest,int * h_edges, pool P,FILE* fp1){
    // int *exclude_S=(int *)malloc(P.numVertex*sizeof(int));
    // int *exclude_T=(int *)malloc(P.numVertex*sizeof(int)); 
    int init = P.target; 
    int count; int N; int* Ad; int num; int i; int j; //int KK;
    while(P.source<P.numVertex){
        P.target=P.source+init;
        while(P.target<P.numVertex){
    count=0; 
    //verify that it is connected directly
    N = degree(h_dest,h_edges,P.source);
    Ad=(int*)malloc(N*sizeof(int));
    memcpy(Ad,&(h_dest[h_edges[P.source]]),N*sizeof(int)); 

    //int check=0;
    for(j = 0; j<N;j++){
        if(Ad[j]==P.target){
            P.h_visited1[P.source]=-1;
            P.h_visited[P.target]=-1;
            count++; 
            break;
        }
    }
    free(Ad);
    //loop while count every component
    while(1){ 
    
        //BFS_host(P);
        BFS_host(P);
        num = *P.h_returned;

        if(num==-1) {
            break;
            }
        count++;

        i=num;
        j=num; 
        P.h_visited[i]=-1;
        P.h_visited1[i]=-1;

        int tempi=0;
        int tempj=0;
        int error =0;
        while((i>-1)||(j>-1)){
            if(i>-1){
                if(i!=num) {
                P.h_visited[i]=-1;
                P.h_visited1[i]=-1;
                }
                tempi=i;
                i = path(P.h_label1,i,h_dest,h_edges,0,P.numVertex);
                if(tempi==i){
                    error=1;
                    break;
                }
            }
            if(j>-1){
                if(j!=num) {
                P.h_visited[j]=-1;
                P.h_visited1[j]=-1;
                }
                tempj=j;
                j = path(P.h_label,j,h_dest,h_edges,1,P.numVertex);
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
    //fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    
    P.target+=4;
    //P.target++;
    for (i=0;i<2*P.numVertex;i++) {
        P.h_visited[i] =0;
        }
        }
    P.source++;
    }

    // free(exclude_S);
    // free(exclude_T);
}

__global__ void BFS_Bqueue(int* p_frontier, int* p_frontier_tail, int* c_frontier, int* c_frontier_tail, int* edges, int* dest, int* label, int* visited, int* label1, int* visited1){
    __shared__ int c_frontier_s[BLOCK_QUEUE_SIZE];
    __shared__ int c_frontier_tail_s, our_c_frontier_tail; 
    int * thisvisited;
    int * thislabel;
    if (threadIdx.x == 0)
        c_frontier_tail_s = 0;
    __syncthreads(); 

    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if (tid < *p_frontier_tail) {
        const int my_vertex = p_frontier[tid]>>1;
        const int my_state = p_frontier[tid]&1;
        if(my_state==1){
            thislabel = label1;
            thisvisited = visited1;
        }
        else{
            thislabel = label;
            thisvisited = visited;
        }
        for (int i = edges[my_vertex]; i < edges[my_vertex+1]; i++) {
            int was_visited=2;
            if(thisvisited[dest[i]]!=-1)
                was_visited = atomicExch(&(thisvisited[dest[i]]), 1);
            if (!was_visited) {
                thislabel[dest[i]] = thislabel[my_vertex] + 1;
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
    __syncthreads();
    if (tid == 0) {
        *p_frontier_tail = atomicExch(c_frontier_tail, 0);
    }
} 

void compute_test(int* h_dest,int * h_edges, pool P){
    int count; int N; int* Ad; int num; int i; int j;

    count=0; 
    //verify that it is connected directly
    N = degree(h_dest,h_edges,P.source);
    Ad=(int*)malloc(N*sizeof(int));
    memcpy(Ad,&(h_dest[h_edges[P.source]]),N*sizeof(int)); 

    //int check=0;
    for(j = 0; j<N;j++){
        if(Ad[j]==P.target){
            P.h_visited[P.source]=-1;
            P.h_visited1[P.target]=-1;
            count++; 
        }
    }
    free(Ad);
    //loop while count every component
    while(1){ 
    
        //BFS_host(P);
        BFS_host(P);
        num = *P.h_returned;
        if(num==-1) {
            break;
            }
        count++;

        i=num;
        j=num; 
        P.h_visited[i]=-1;
        P.h_visited1[i]=-1;


        while((i>-1)||(j>-1)){
            if(i>-1){
                if(i!=num) {
                P.h_visited[i]=-1;
                P.h_visited1[i]=-1;
                }
                // tempi=i;
                i = path(P.h_label1,i,h_dest,h_edges,0,P.numVertex);
                // if(tempi==i){
                //     error=1;
                //     break;
                // }
            }
            if(j>-1){
                if(j!=num) {
                P.h_visited[j]=-1;
                P.h_visited[j]=-1;
                }
                // tempj=j;
                j = path(P.h_label,j,h_dest,h_edges,1,P.numVertex);
                // if(tempj==j){
                //     error=1;
                //     break;
                // }
            }
        }
        // if(error==1){
        //     count--;
        //     break;
        // }
    }
    //record count
    //fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    //printf("[%d, %d] %d\n", P.source, P.target, count);
    
}
