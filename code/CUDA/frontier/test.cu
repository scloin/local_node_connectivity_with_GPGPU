#include <revise4t.hpp>
#define BLOCK_SIZE 512
#define BLOCK_QUEUE_SIZE 512

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
    int *h_data;
    h_data = (int*)malloc((elen+dlen)*sizeof(int)); 
    //CUDA_CHECK(cudaMallocHost((void**)&h_data, (elen+dlen)*sizeof(int)));
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
    int* devmem; 
    int* d_edges;
    CUDA_CHECK(cudaMalloc((void**)&d_edges, (elen+dlen)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem, (6*(elen-1)+2)*sizeof(int)));

    pool P0 = init_pool(elen, dlen, devmem, d_edges);
    CUDA_CHECK(cudaMemcpy(d_edges, h_edges, (dlen+elen)*sizeof(int), cudaMemcpyHostToDevice)); 

    ///////////////////////////////////////////////////////
    /*compute*/
    P0.target=P0.source+1;
    compute_test(h_dest,h_edges,P0);



    ///////////////////////////////////////////////////////
    /*free*/
    CUDA_CHECK(cudaFree(devmem));
    CUDA_CHECK(cudaFree(d_edges)); 
    CUDA_CHECK( cudaStreamDestroy(P0.stream));
    free(h_data);
    CUDA_CHECK(cudaFreeHost(P0.h_label));
    free(P0.h_returned);

    cudaDeviceReset();
    return 0;
} 