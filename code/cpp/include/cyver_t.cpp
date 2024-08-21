#include <cyver_t.hpp>



int path(int* exclude, int num, int* h_dest, int* h_edges, int state, int ENUM){
    int numT=num;
    int check=0;
    int i=0;

    int N = degree(h_dest,h_edges,num);
    
    if(exclude[numT]==1){
        if(state==0) {//printf("next is target ");
            return -2;
        }else {//printf("next is source ");
            return -1;
        }//break;
    }
    else{
    int* AA=(int*)malloc(N*sizeof(int));
    memcpy(AA,&(h_dest[h_edges[num]]),N*sizeof(int));
        for (i = 0; i < ENUM; i++){
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
    if (check!=1)  i=numT;

    return i;
}

int degree(int* dest,int* edges,int source){
    
    return edges[source+1]-edges[source];
}

void _bidirectional_pred_succ(int* dest,int* edges,int source,int target,int* label_S,int* label_T,int* returned, int ENUM){

    // does BFS from both source and target and meets in the middle
    // excludes nodes in the container "exclude" from the search

    if (target == source) {
        //printf("source and target can't be the same");
        *returned = -1;
    }

    int frontier[4][int(ENUM)];
    
    int *c_frontier_S = frontier[0];
    int *p_frontier_S = frontier[1];
    int *c_frontier_T = frontier[2];
    int *p_frontier_T = frontier[3];  

    int c_frontier_S_tail = 0;
    int p_frontier_S_tail = 0;
    int c_frontier_T_tail = 0;
    int p_frontier_T_tail = 0;  

    insert_frontier(source, p_frontier_S, &p_frontier_S_tail);
    insert_frontier(target, p_frontier_T, &p_frontier_T_tail);

    int level = 0;
    int check = 0;

    label_S[source] = 0;
    label_T[target] = 0;

    //int *tmp;
    int num =-1;
    
    while ((p_frontier_S_tail > 0)&&(check==0)) {
        level+=1;
    if ((level%2)==0){
        for (int f = 0; f < p_frontier_S_tail; f++) {
            int c_vertex = p_frontier_S[f];
            for (int i = edges[c_vertex]; i < edges[c_vertex+1]; i++) {
                if (label_S[dest[i]] == -1) {
                    insert_frontier(dest[i], c_frontier_S, &c_frontier_S_tail);
                    label_S[dest[i]] = label_S[c_vertex] + 1;
                }
                if (label_T[dest[i]] > -1) {
                    //printf("\n\n====CHECK!(%d)====\n", dest[i]);
                    check=1;
                    num=dest[i];
                    //return num;
                }
                if (check==1) break;
            }
            if (check==1) break;
        }
        if (check==1) continue;
        int *tmp = c_frontier_S;
        c_frontier_S = p_frontier_S;
        p_frontier_S = tmp;

        p_frontier_S_tail = c_frontier_S_tail;
        c_frontier_S_tail = 0;
        }
    else{
        for (int f = 0; f < p_frontier_T_tail; f++) {
            int c_vertex = p_frontier_T[f];
            for (int i = edges[c_vertex]; i < edges[c_vertex+1]; i++) {
                //printf("%d", dest[i]);
                if (label_T[dest[i]] == -1) {
                    insert_frontier(dest[i], c_frontier_T, &c_frontier_T_tail);
                    label_T[dest[i]] = label_T[c_vertex] + 1;
                }
                if (label_S[dest[i]] > -1) {
                    //printf("\n\n====CHECK!(%d)====\n", dest[i]);
                    check=1;
                    num= dest[i];
                    //return num;
                }
                if (check==1) break;
            }
            if (check==1) break;
        }
        if (check==1) continue;
        int *tmp = c_frontier_T;
        c_frontier_T = p_frontier_T;
        p_frontier_T = tmp;

        p_frontier_T_tail = c_frontier_T_tail;
        c_frontier_T_tail = 0; 
    }

    if (check==1) break;
    }
    //free(p_frontier_S);
    //free(p_frontier_T);
    //free(c_frontier_S);
    //free(c_frontier_T);
    //free(tmp);
    *returned = num;
}

void insert_frontier(int source, int* frontier, int* frontier_tail)
{
    frontier[*frontier_tail] = source;
    (*frontier_tail)++;
}

void compute_test(int* h_dest,int * h_edges, pool P){
    int *exclude_S=(int *)malloc(P.numVertex*sizeof(int));
    int *exclude_T=(int *)malloc(P.numVertex*sizeof(int)); 
    int count; int N; int* Ad; int num; int i; int j; int KK;
    
    count=0; 
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
        }
    }
    free(Ad);
    //loop while count every component
    while(1){ 
    
        //BFS_host(P);
        //_bidirectional_pred_succ(P.dest, P.edges);
        num = *P.h_returned;
        if(num==-1) {
            break;
            }
        count++;

        i=num;
        j=num; 
        P.h_visited[i<<1]=-1;
        P.h_visited[(i<<1)+1]=-1;

        for (KK=0; KK<P.numVertex;KK++){
                exclude_S[KK] = P.h_label[KK*2];
                exclude_T[KK] = P.h_label[KK*2+1];
            } 

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
    //fprintf(fp1,"[%d, %d] %d\n", P.source, P.target, count);
    
    free(exclude_S);
    free(exclude_T);
}