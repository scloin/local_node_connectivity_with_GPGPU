#include <revise5t.hpp>
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
    int* devmem1; 
    int* devmem2;
    int* devmem3;
    int* d_edges;
    CUDA_CHECK(cudaMalloc((void**)&d_edges, (elen+dlen)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem, (6*(elen-1)+2)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem1, (6*(elen-1)+2)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem2, (6*(elen-1)+2)*sizeof(int)));
    CUDA_CHECK(cudaMalloc((void**)&devmem3, (6*(elen-1)+2)*sizeof(int)));

    pool P0 = init_pool(elen, dlen, devmem, d_edges);
    pool P1 = init_pool(elen, dlen, devmem1, d_edges);
    pool P2 = init_pool(elen, dlen, devmem2, d_edges);
    pool P3 = init_pool(elen, dlen, devmem3, d_edges);
    CUDA_CHECK(cudaMemcpy(d_edges, h_edges, (dlen+elen)*sizeof(int), cudaMemcpyHostToDevice)); 
    thread t0, t1, t2, t3;
    FILE* fp1 = fopen("result/no_frontier.txt","w"); 

    ///////////////////////////////////////////////////////
    /*compute*/
    P0.target=P0.source+1;
    P1.target=P0.source+2;
    P2.target=P0.source+3;
    P3.target=P0.source+4;
    P1.source=P0.source; 
    P2.source=P0.source;
    t0=thread{compute,h_dest,h_edges,P0,fp1};
    t1=thread{compute,h_dest,h_edges,P1,fp1};
    t2=thread{compute,h_dest,h_edges,P2,fp1};
    t3=thread{compute,h_dest,h_edges,P3,fp1};
    //compute(h_dest,h_edges,P0,fp1);
    t0.join();
    t1.join();
    t2.join();
    t3.join();

    fclose(fp1);

    ///////////////////////////////////////////////////////
    /*free*/
    CUDA_CHECK(cudaFree(devmem));
    CUDA_CHECK(cudaFree(devmem1));
    CUDA_CHECK(cudaFree(devmem2));
    CUDA_CHECK(cudaFree(devmem3));
    CUDA_CHECK(cudaFree(d_edges)); 
    CUDA_CHECK( cudaStreamDestroy(P0.stream));
    CUDA_CHECK( cudaStreamDestroy(P1.stream));
    CUDA_CHECK( cudaStreamDestroy(P2.stream));
    CUDA_CHECK( cudaStreamDestroy(P3.stream));
    free(h_data);
    free(P0.h_label);
    free(P0.h_returned);
    free(P1.h_label);
    free(P1.h_returned);
    free(P2.h_label);
    free(P2.h_returned);
    free(P3.h_label);
    free(P3.h_returned);
    //cudaFreeHost(h_data);
    cudaDeviceReset();
    return 0;
} 