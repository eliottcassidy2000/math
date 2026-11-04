#include <iostream>
#include <cstdint>
using namespace std;
int main(){for(auto tf:{pair<int,uint64_t>{1,2710066},{1,3035358},{1,3186353},{1,3568813},{2,1355033},{2,1517679}}){int t=tf.first;uint64_t f=tf.second,cut=f<<t;int n=t+(64-__builtin_clzll(f));uint64_t all=(1ULL<<n)-1,h=0;for(uint64_t before=0;before<=all;before++){uint64_t after=all^before;if(before && !(cut>>(63-__builtin_clzll(before))&1))continue;if(after && (cut>>__builtin_ctzll(after)&1))continue;h++;}if(h!=4611841)return 1;cout<<t<<' '<<f<<' '<<n+1<<' '<<h<<'\n';}}
