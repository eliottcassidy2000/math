#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
int main(){
 constexpr int n=21;
 constexpr uint32_t L=1<<n,V=L-1;
 const int S[][2]={{0,1},{1,0},{1,1},{1,2},{2,0},{2,1},{2,2},{3,0},{4,1},{4,2}};
 std::array<uint32_t,n> A{};
 for(int u=0;u<n;u++){
   int a=u%7,b=u/7;
   for(auto& s:S) A[u] |= 1<<(((a+(1<<b)*s[0])%7)+7*((b+s[1])%3));
 }
 // Avoid sets containing identity: the DP computes suffixes after the fixed first vertex.
 const int m=n-1; const uint32_t W=1<<m;
 std::vector<uint64_t> B((size_t)m*W);
 for(int u=0;u<m;u++) B[(size_t)(1<<u)*m+u]=1;
 for(uint32_t x=1;x<W;x++){
   if ((x&(x-1))==0) continue;
   for(uint32_t ubit=x;ubit;ubit&=ubit-1){
     int u=std::countr_zero(ubit);
     const uint32_t y=x^(1u<<u);
     uint64_t ans=0;
     for(uint32_t vbit=y&(A[u+1]>>1);vbit;vbit&=vbit-1) ans+=B[(size_t)y*m+std::countr_zero(vbit)];
     B[(size_t)x*m+u]=ans;
   }
 }
 uint64_t total=0;
 for(int u=1;u<n;u++)if((A[0]>>u)&1){
   uint64_t ans=B[(size_t)(W-1)*m+u-1]; total+=ans;
   std::cout << u%7 << ","<< u/7 << ":" << ans <<"\n";
 }
 std::cout <<"rooted "<<total<<" total "<<21*total<<"\n";
}
