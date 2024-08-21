#include <cyver_t.hpp>





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
    int ENUM=elen-1;
    printf("N : %d\n", elen-1);
 
    int *exclude_S = (int*)malloc((ENUM)*sizeof(int));
    int *exclude_T = (int*)malloc((ENUM)*sizeof(int));
    int source=0;
    int target=0;

    FILE* fp1 = fopen("./result/casecpp.txt","w"); //test파일을 w(쓰기) 모드로 열기
    std::list<int> list1;
    while(source<(elen-1)){

    target=source+1;

    while(target<(elen-1)){  
        for (int i = 0; i < (elen-1); i++){
            exclude_S[i] = -1;
            exclude_T[i] = -1;
            }
    
    int num=-1;
    _bidirectional_pred_succ(h_dest,h_edges,source,target,exclude_S,exclude_T,&num, ENUM);
    if(num==-1) {
        //printf("%d->%d; none\n", source,target);
        fprintf(fp1,"%d->%d; none\n", source,target);
        target++;
        continue;
    }
    else if(num==source || num==target) {
        //printf("%d %d \n", source,target);
        fprintf(fp1,"%d %d \n", source,target);
        target++;
        continue;
    }
    list1.push_front(source);
    list1.push_back(target);
    std::list<int>::iterator begin_iter = list1.begin();
    std::list<int>::iterator end_iter = list1.end();
    end_iter--; begin_iter++;
    list1.insert(end_iter, num);
    begin_iter--;
    //int *AA;
    int i=num;
    while((i>-1)){
        if(i!=num) {
        list1.insert(end_iter, i);}
        i = path(exclude_T,i,h_dest,h_edges,0,ENUM);
        }

    i=num;
    while((i>-1)){
        if(i!=num) {
        list1.insert(begin_iter, i);
        begin_iter--;}
        i = path(exclude_S,i,h_dest,h_edges,1,ENUM);
        }

    while (list1.empty()==0) {
        i=list1.front();
        fprintf(fp1,"%d ",i); //문자열 입력
        //printf("%d ",i);
        list1.pop_front();

    }

    fprintf(fp1,"\n");
    //printf("\n");
    list1.clear();
    target++;
    }
    source++;}
    //fclose(fp1);
    free(h_dest);
    free(h_edges);
    //free(h_label);
    free(exclude_T);
    free(exclude_S);
    return 0;
}