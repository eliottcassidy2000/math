#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
uint64_t calc(int k,uint64_t A,uint64_t b){int n=2*k,N=1<<k,adj[24]={},z=0;auto edge=[&](int u,int v,int yes){if(yes)adj[u]|=1<<v;else adj[v]|=1<<u;};for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){int a=A>>z&1,c=b>>z&1;edge(2*i,2*j,a);edge(2*j+1,2*i+1,a);edge(2*i,2*j+1,c);edge(2*j,2*i+1,c);}static uint64_t dp[4096][24];fill(&dp[0][0],&dp[N-1][24],0);for(int i=0;i<n;i++)dp[1<<(i/2)][i]=1;for(int mask=1;mask<N;mask++){int av=0;for(int j=0;j<k;j++)if(!(mask>>j&1))av|=3<<(2*j);for(int i=0;i<n;i++){auto v=dp[mask][i];if(v){int next=adj[i]&av;while(next){int j=__builtin_ctz(next);next&=next-1;dp[mask|(1<<(j/2))][j]+=v;}}}}uint64_t h=0;for(int i=0;i<n;i+=2)h+=dp[N-1][i];return h;}
int main(int argc,char**argv){if(argc!=4)return 2;int r=stoi(argv[1]),k=r+1;if(k>10)return 2;uint64_t A=stoull(argv[2]),b=stoull(argv[3]),a0=0,b0=0;vector<int> free;int q=0,z=0;for(int i=0;i<k;i++)for(int j=i+1;j<k;j++,z++){if(j==r)free.push_back(z);else{a0|=((A>>q)&1)<<z;b0|=((b>>q)&1)<<z;q++;}}uint64_t best=0,ba=0,bb=0;int count=0;for(int c=0;c<(1<<(2*r));c++){uint64_t a=a0,d=b0;for(int i=0;i<r;i++){a|=(uint64_t)(c>>i&1)<<free[i];d|=(uint64_t)(c>>(r+i)&1)<<free[i];}auto h=calc(k,a,d);if(h>best){best=h;ba=a;bb=d;count=0;}if(h==best)count++;}cout<<"{\"k\":"<<k<<",\"A\":"<<ba<<",\"b\":"<<bb<<",\"h\":"<<best<<",\"max_extension_count\":"<<count<<",\"extensions_checked\":"<<(1<<(2*r))<<"}\n";}
