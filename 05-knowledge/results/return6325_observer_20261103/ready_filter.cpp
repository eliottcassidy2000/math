#include <iostream>
#include <vector>
#include <cstdint>
int main(){
 int n,b,index;uint64_t target;size_t processed=0;
 while(std::cin>>n>>b>>index>>target){
  std::vector<unsigned> T(n);for(auto &x:T)std::cin>>x;
  int size=1<<n,full=size-1;std::vector<uint64_t> dp(size*n),hs(size);
  for(int v=0;v<n;v++)dp[(1<<v)*n+v]=1;
  for(int S=1;S<size;S++)for(int v=0;v<n;v++){
   auto z=dp[S*n+v];hs[S]+=z;
   if(!z)continue;unsigned rest=T[v]&~S;
   while(rest){unsigned bit=rest&-rest;rest-=bit;int j=__builtin_ctz(bit);dp[(S|bit)*n+j]+=z;}
  }
  hs[0]=1;uint64_t F=0;for(int S=1;S<size;S+=2)F+=dp[S*n]*hs[full^S];
  if(hs[full]+F-b==target)std::cout<<index<<" "<<hs[full]<<" "<<hs[full^1]<<" "<<F<<"\n";
  processed++;
 }
 std::cerr<<"processed "<<processed<<"\n";
}
