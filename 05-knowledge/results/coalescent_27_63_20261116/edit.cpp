#include <array>
#include <functional>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
int main(){
 const int n=9,N=1<<n;
 std::array<int,n> a{};
 std::vector<std::pair<int,int>> edges;
 for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) {
  edges.emplace_back(i,j);
  if(i/3==j/3 && i%3==0 && j%3==2) a[j]|=1<<i;
  else a[i]|=1<<j;
 }
 auto count=[&](){
  int dp[N][n]={};
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(dp[s][i]){
   for(int b=a[i]&~s;b;b&=b-1){int j=__builtin_ctz(b);dp[s|1<<j][j]+=dp[s][i];}
  }
  return std::accumulate(dp[N-1],dp[N-1]+n,0);
 };
 if(count()!=27)return 1;
 std::vector<std::pair<int,int>> flips;
 for(int k=0;k<=5;k++){
  unsigned long long hit=0;
  auto rec=[&](auto&&self,int left,int p)->void{
   if(left==0){if(count()==63){
     if(hit==0){for(auto [u,v]:flips)std::cout<<u<<"-"<<v<<" ";std::cout<<" rows";for(auto r:a)std::cout<<" "<<r;std::cout<<"\n";}
     ++hit;
   }return;}
   for(int i=p;i<=int(edges.size())-left;i++){
    auto [u,v]=edges[i];a[u]^=1<<v;a[v]^=1<<u;
    flips.emplace_back(u,v);self(self,left-1,i+1);flips.pop_back();
    a[u]^=1<<v;a[v]^=1<<u;
   }
  };
  rec(rec,k,0);
  std::cout<<k<<" flips: "<<hit<<" hits\n"<<std::flush;
  if(hit)break;
 }
}
