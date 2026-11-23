#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
// Subset path DP, compressed by deleting the last-vertex bit from each column.
constexpr int N=23, M=1<<N;
using u64=uint64_t;
using A=std::array<unsigned,N>;
int prod(int x,int y){int b=x/7;return (x%7+(1<<b)*(y%7))%7+7*((b+y/7)%3);}
A make(bool dial,bool swapped){
  A R{};
  for (int i=0;i<21;i++)for(int j:{1,2,3,14,8,9,10,11,12,13}){
      int k=prod(i,j);
      if(dial)R[k]|=1U<<i;else R[i]|=1U<<k;
  }
  for(int i=0;i<2;i++){
    int a=i^(int)swapped;int v=21+i;
    for(int u=0;u<21;u++)if(u==a)R[v]|=1U<<u;else R[u]|=1U<<v;
  }
  R[21]|=1U<<22;
  return R;
}
u64 height(A R){
  A I{};for(int i=0;i<N;i++)for(int j=0;j<N;j++)if(R[i]>>j&1)I[j]|=1U<<i;
  std::vector<u64> paths(size_t(N)*(M>>1));std::array<u64*,N> p{};
  for(int v=0;v<N;v++){p[v]=paths.data()+size_t(v)*(M>>1);p[v][0]=1;}
  auto compressed=[](unsigned s,int j) {unsigned l=(1U<<j)-1;return (s&l)|((s>>1)&~l);};
  for(unsigned s=1;s<M;s++){
    if(!(s&(s-1)))continue;
    for(unsigned u=s;u;u&=u-1){
      int v=__builtin_ctz(u);unsigned t=s^(1U<<v);
      u64 sum=0;
      for(unsigned r=t&I[v];r;r&=r-1){int j=__builtin_ctz(r);sum+=p[j][compressed(t,j)];}
      p[v][compressed(t,v)]=sum;
    }
  }
  u64 h=0;for(auto x:p)h+=x[(M>>1)-1];
  return h;
}
int main(){
  for (bool d:{false,true})
    for(bool s:{false,true})std::cout<<d<<" "<<s<<" "<<height(make(d,s))<<std::endl;
}
