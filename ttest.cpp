#include<iostream>
#include<thread>
#include<stdio.h>
#include <vector>

using namespace std;
 
void worker(vector<int>::iterator start, vector<int>::iterator end,int* result) {
    int sum = 0;
    for (auto itr = start; itr < end; ++itr) sum += *itr;
    *result = sum;
    // 쓰레드의 id 를 구한다.
    thread::id this_id = std::this_thread::get_id();
    printf("쓰레드 %p 에서 %d 부터 %d 까지 계산한 결과 : %d \n", this_id, *start,*(end - 1), sum);
    //cout을 사용하면 << 시 스레가 바뀌면서 메세지가 뒤바껴서 나올수 있음 
}
 
int main()
{
    int ans = 0;
    vector<int> v;
    vector<int> sum(5);
    vector<thread> t;
    thread t1;
    for (int i = 0; i < 10000; i++) v.push_back(i);
    for (int i = 0; i < 4; i++) {
        t.push_back(thread{worker,v.begin() + i * 2500,v.begin() + (i + 1) * 2500,&sum[i] });
    }// 첫 번째 스레드는 0~2499 두번째는 2500~4999 ... callable의 개념으로 호출
    int i=3;
    t1=thread{worker,v.begin() + i * 2500,v.begin() + (i + 1) * 2500, &sum[4]};
    t1.join(); 
    printf("뭐임? %d\n", sum[i]);
    for (int i = 0; i < 4; i++) t[i].join();
    for (int i = 0; i < 4; i++) ans += sum[i];
    cout << "전체 합 : " << ans;
    return 0;
}