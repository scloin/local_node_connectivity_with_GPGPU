//nvcc -o test test/cu -I.. -arch=sm_35 -rdc=true
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <list>
#include <common/common.h>
#include <common/common_string.h> 

#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 32 

__global__ void add(int* A, int* B, int* C);
//__global__ void filter(int* A, int* C);
//__global__ void mainstream(int* A, int* B, int* C);

__global__ void BFSmain(int* in,int* out, int* edges, int*dest, int N,int* kk);
__global__ void matmul(int* input, int* output, int* edges, int*dest, int NUM, int T);
__global__ void filter(int* input, int* output, int NUM);
__global__ void checking(int* input, int* input1, int* kk, int NUM);

int path(int* exclude, int num, int* h_dest, int* h_edges, int state,int elen);
int degree(int* dest,int* edges,int source);

int main(){
    //init&alloc
        // int A[]={0,1,2,3,4,5,6,7};
        // int B[]={0,1,2,3,4,5,6,7};
        // //int C[]={0,0,0,0,0,0,0,0};
        // int* C=(int*)calloc(8,sizeof(int));
        // int* Ac, *Bc, *Cc;
        // CUDA_CHECK(cudaMalloc((void**)&Ac, 8*sizeof(int))); 
        // CUDA_CHECK(cudaMalloc((void**)&Bc, 8*sizeof(int)));
        // CUDA_CHECK(cudaMalloc((void**)&Cc, 8*sizeof(int)));
        
        // CUDA_CHECK(cudaMemcpy(Ac, A, 8*sizeof(int), cudaMemcpyHostToDevice)); 
        // CUDA_CHECK(cudaMemcpy(Bc, B, 8*sizeof(int), cudaMemcpyHostToDevice)); 
        // CUDA_CHECK(cudaMemcpy(Cc, C, 8*sizeof(int), cudaMemcpyHostToDevice)); 

    //read file&gengraph
        FILE* fp = fopen("./graph/graph.txt","r");
        int dlen; int elen; 
        int * h_dest; int * h_edges; 
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
        int *d_edges, *d_dest;
        CUDA_CHECK(cudaMalloc((void**)&d_edges, elen*sizeof(int)));
        CUDA_CHECK(cudaMalloc((void**)&d_dest, dlen*sizeof(int)));
        CUDA_CHECK(cudaMemcpy(d_edges, h_edges, elen*sizeof(int), cudaMemcpyHostToDevice)); 
        CUDA_CHECK(cudaMemcpy(d_dest, h_dest, dlen*sizeof(int), cudaMemcpyHostToDevice)); 
    int num=elen-1;
    int* h_in, *d_in, *h_out, *d_out;
    int* kk;
    int source; int target; int count;
    h_in=(int*)calloc(2*num,sizeof(int));
    h_out=(int*)calloc(2*num,sizeof(int));
    CUDA_CHECK(cudaMalloc((void**)&d_in, 2*num*sizeof(int))); 
    CUDA_CHECK(cudaMalloc((void**)&d_out, 2*num*sizeof(int))); 
    CUDA_CHECK(cudaMalloc((void**)&kk, 2*sizeof(int)));

    FILE* fp1 = fopen("result/test4.txt","w");

    source=0; target=12; 
    while(source<num){
    target=source+1;
    while(target<num){
    count=0; int meet = -1;

    int exclude_S[num];
    int exclude_T[num];
    do{
    h_in[source]=1; h_in[target+num]=1; h_out[source]=-1; h_out[target+num]=-1;
    CUDA_CHECK(cudaMemset(kk, 0, 2*sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_in, h_in, 2*num*sizeof(int), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(d_out, h_out, 2*num*sizeof(int), cudaMemcpyHostToDevice)); 

    //while()
    BFSmain<<<1,2>>>(d_in, d_out, d_edges, d_dest, num,kk);
    CUDA_CHECK(cudaMemcpy(h_out, d_out, 2*num*sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_in, d_in, 2*num*sizeof(int), cudaMemcpyDeviceToHost));
    int min = num;
    meet = -1; 

    for(int k=0;k<num;k++){
        if(h_in[k]==1&&h_in[k+num]==1){
            if((min>h_out[k]+h_out[k+num])&&(h_out[k]+h_out[k+num]>=0)){

                if((k==source || k==target)&&(h_out[k]+h_out[k+num]==0)){
                    count+=1;
                    h_in[source]=-1;
                    h_in[source+num]=-1;
                    h_in[target]=-1;
                    h_in[target+num]=-1;
                    continue;
                }
                else if((k==source || k==target)){
                    h_in[source]=-1;
                    h_in[source+num]=-1;
                    h_in[target]=-1;
                    h_in[target+num]=-1;
                    continue;
                }
                min = h_out[k]+h_out[k+num];
                meet =k;
            }
        }}
    for(int k=0;k<num;k++){
        if(h_out[k]==0) h_out[k]=-1;
        else if(h_out[k]==-1) h_out[k]=0;
        if(h_out[k+num]==0) h_out[k+num]=-1;
        else if(h_out[k+num]==-1) h_out[k+num]=0;
    }
    
    if(meet==-1){break;}
    else{
    count++;
    std::list<int> list1;    
    std::list<int>::iterator begin_iter = list1.begin();
    std::list<int>::iterator end_iter = list1.end();
    list1.insert(end_iter, meet);
    
    int i=meet; 
    //printf("%d\n",i);
    for (int KK=0; KK<num;KK++){
            exclude_S[KK] = h_out[KK];
            exclude_T[KK] = h_out[KK+num];
        } 

    while((i>-1)){
        if(i!=meet) {
        list1.insert(end_iter, i);}
        i = path(exclude_T,i,h_dest,h_edges,0,num);
        } 

    i=meet;
    while((i>-1)){
        if(i!=meet) {
        list1.insert(begin_iter, i);
        begin_iter--;}
        i = path(exclude_S,i,h_dest,h_edges,1,num); 

        }
    
    while (list1.empty()==0) {
        i=list1.front();
        h_in[i]=-1;
        h_in[i+num]=-1;
        list1.pop_front(); 
        //if(source==10&&target==12){printf("%d ",i);}
    } 

    list1.clear();
    //if(source==10&&target==12) printf(" [%d]\n",count);
    }
    for (int i=0;i<2*num;i++) {
        h_out[i]=0;
        if(h_in[i]!=-1) h_in[i] =0;
    }
    }while(meet!=-1);
    fprintf(fp1,"[%d, %d] %d\n", source, target, count);
    target++;
        for (int i=0;i<2*num;i++) {
        h_out[i]=0;
        h_in[i]=0;
        }
    } 
    source++;
    }
    //printf("%d\n", count);
    fclose(fp1);

    free(h_in); free(h_out); free(h_dest); free(h_edges);
    CUDA_CHECK(cudaFree(d_in));
    CUDA_CHECK(cudaFree(d_out));
    CUDA_CHECK(cudaFree(d_dest));
    CUDA_CHECK(cudaFree(d_edges));
    CUDA_CHECK(cudaFree(kk));
    return 0;
}

/*
    __global__ void add(int* A,int *B,int *C){
        const int tid = blockIdx.x*blockDim.x + threadIdx.x;
        
        if(tid<8){
            C[tid]=A[tid]*B[tid];
        }__syncthreads();
    }

    __global__ void filter(int* A,int *C){
        const int tid = blockIdx.x*blockDim.x + threadIdx.x;
        if(tid<8){
            C[tid]=C[tid]%9;
            A[tid]=C[tid];
        }__syncthreads();
    }

    __global__ void mainstream(int* A,int *B,int *C){
        int i=0;
        while(i<5){
        add<<<1,8>>>(A,B,C);
        //__syncthreads();
        cudaDeviceSynchronize();
        filter<<<8,1>>>(A,C);
        cudaDeviceSynchronize();
        i++;
        for(int i=0; i<8;i++) printf("%d ",C[i]);
        printf("\n");
        }
    }
*/
__global__ void BFSmain(int* in,int* out, int* edges, int*dest, int N, int* kk){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int ccheck=1; int check=0;
    int num_blocks = ((N*N)+BLOCK_SIZE-1) / BLOCK_SIZE;
    int num_blocks1 = (N+BLOCK_SIZE-1) / BLOCK_SIZE;
    int* in1=&in[N]; int* out1=&out[N];
    while(kk[0]==0){
        check=kk[1];
        kk[1]=0;
        if(tid==0){
            matmul<<<num_blocks,BLOCK_SIZE>>>(in,out, edges, dest, N, ccheck);
            filter<<<num_blocks1,BLOCK_SIZE>>>(in,out, N);}
        else if(tid==1){
            matmul<<<num_blocks,BLOCK_SIZE>>>(in1,out1, edges, dest, N, ccheck);
            filter<<<num_blocks1,BLOCK_SIZE>>>(in1,out1, N);}
        cudaDeviceSynchronize();

        //filter<<<num_blocks1,BLOCK_SIZE>>>(in,out, N);
        //filter<<<num_blocks1,BLOCK_SIZE>>>(in1,out1, N);
        //cudaDeviceSynchronize();
        if(tid==0)
        checking<<<num_blocks1,BLOCK_SIZE>>>(in,in1,kk,N);
        cudaDeviceSynchronize();
        ccheck++;
        if(check==kk[1]) kk[0]=-1;
    }
}

__global__ void matmul(int* input, int* output, int* edges, int*dest, int NUM,int T){
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int N = tid/NUM; int M = tid%NUM;
    int k=0;
    if(M<(edges[N+1]-edges[N])){
        k=input[dest[edges[N]+M]];
        if(output[N]==0&&k!=0&&input[N]!=-1&&k!=-1){
            atomicExch(&output[N],T);
        }
    }
    __syncthreads();
}

__global__ void filter(int* input, int* output, int NUM){
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    //int check0=0;
    //int dd[1];
    //dd[0]=1;
    if(tid<NUM){
        //printf("              wht?%d[%d]\n",tid,output[tid]);
        if(output[tid]==0) {
            //printf("                  if %d\n",tid);
            //atomicAdd(&check0,1);
        }
        else if(input[tid]==0){
            input[tid]=1;
            //printf("                  elseif %d\n",tid);
        }
        else if(input[tid]==1){
            //input[tid]=0;
            //printf("                  else %d\n",tid);
            }
    }
    // if(tid==1){
    //     for(int i=0; i<NUM;i++) printf("%d ",input[i]);
    //         printf("\n");}
    __syncthreads();
}

__global__ void checking(int* input, int* input1, int* kk, int NUM){
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if(tid<NUM){
        if(input[tid]==input1[tid]&&input1[tid]==1){
            kk[0]=tid;
        }
        else if(input[tid]==0||input1[tid]==0){
            atomicAdd(&kk[1],input[tid]+input1[tid]);

        }
    }
    __syncthreads();
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
