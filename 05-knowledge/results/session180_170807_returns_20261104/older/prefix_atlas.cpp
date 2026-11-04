#include <array>
#include <iostream>
#include <map>
#include <tuple>
#include <cstdint>
using namespace std;
int main(){map<pair<int,int>,pair<int,uint32_t>> atlas;long long count=0;
for(int n=1;n<=7;n++){int N=1<<n,m=n*(n-1)/2; for(uint32_t code=0;code<(1u<<m);code++){int incoming[7]={};int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(code>>k&1)incoming[j]|=1<<i;else incoming[i]|=1<<j;
int dp[128][7]={},h[128]={};h[0]=1;for(int s=1;s<N;s++){int ends=s;while(ends){int v=__builtin_ctz(ends);ends&=ends-1;int t=s^(1<<v);if(!t)dp[s][v]=1;else{int p=t&incoming[v];while(p){int u=__builtin_ctz(p);p&=p-1;dp[s][v]+=dp[t][u];}}h[s]+=dp[s][v];}}
int p2=0;for(int s=0;s<N;s++)p2+=h[s]*h[(N-1)^s];if(p2<=198)atlas.emplace(make_pair(h[N-1],p2),make_pair(n,code));count++;}
cerr<<n<<" complete, pairs "<<atlas.size()<<"\n";}
cout<<"{\"labelled_rows\":"<<count<<",\"max_p2\":198,\"pairs\":[";bool first=true;for(auto [p,w]:atlas){if(!first)cout<<",";first=false;cout<<"{\"h\":"<<p.first<<",\"p2\":"<<p.second<<",\"n\":"<<w.first<<",\"code\":"<<w.second<<"}";}cout<<"]}\n";}
