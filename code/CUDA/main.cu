#include <revise4t.hpp>
#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 32 

using namespace std;

int main(){

    ///////////////////////////////////////////////////////
    /*read file*/
    FILE* fp = fopen("./graph/graph.txt","r");
    int dlen;
    int elen; 
    int * h_dest;
    int * h_edges; 

    fscanf(fp, "%d ", &dlen);
    fscanf(fp, "%d ", &elen); 
    printf("%d\n",elen);
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

    ///////////////////////////////////////////////////////
    /*alloc & init*/
    int* devmem; int* devmem1;
    int* d_edges;
    CUDA_CHECK(cudaMalloc((void**)&d_edges, (elen+dlen)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem, (6*(elen-1)+2)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem1, (6*(elen-1)+2)*sizeof(int)));

    pool P0 = init_pool(elen, dlen, devmem, d_edges);
    pool P1 = init_pool(elen, dlen, devmem1, d_edges);
    CUDA_CHECK(cudaMemcpy(d_edges, h_edges, (dlen+elen)*sizeof(int), cudaMemcpyHostToDevice)); 

    FILE* fp1 = fopen("result/addthread.txt","w"); 

    ///////////////////////////////////////////////////////
    /*compute*/
    while(P0.source<P0.numVertex){
        P0.target=P0.source+1;
        P1.source=P0.source;
        while(P0.target<P0.numVertex){
            P1.target=P0.target+1;
            if(P1.target>=P0.numVertex){
                compute(h_dest,h_edges,P0,fp1);
            }

            else{
                compute(h_dest,h_edges,P0,fp1);
                compute(h_dest,h_edges,P1,fp1);
            }
            P0.target+=2;

            for (int i=0;i<2*P0.numVertex;i++) {
                P0.h_visited[i] =0;
                P1.h_visited[i] =0;
                }
        } 

    P0.source++;}
    fclose(fp1);

    ///////////////////////////////////////////////////////
    /*free*/
    CUDA_CHECK(cudaFree(devmem));
    CUDA_CHECK(cudaFree(devmem1));
    CUDA_CHECK(cudaFree(d_edges)); 
    free(h_data);
    free(P0.h_label);
    free(P0.h_returned);
    free(P1.h_label);
    free(P1.h_returned);
    cudaDeviceReset();
    return 0;
} 