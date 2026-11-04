#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
using U=unsigned;
using I=std::uint64_t;
int main(){
 std::array<std::array<I,6>,10> counts{};
 U n; I code;
 while(std::cin>>n>>code){
  if(n<2||n>9)throw std::runtime_error("bank order");
  U full=(1u<<n)-1,bit=0;std::vector<U>a(n);
  for(U u=0;u<n;++u)for(U v=u+1;v<n;++v,++bit)
   if(code>>bit&1)a[u]|=1u<<v;else a[v]|=1u<<u;
  std::vector<U>g((full+1)*n),s(full+1),d(n),h(n),general(n);
  U boundary=0;
  for(U m=1;m<=full;++m){
   U endmask=0;
   for(U todo=m;todo;todo&=todo-1){
    U u=__builtin_ctz(todo),rest=m^(1u<<u),ends=rest?0:1u<<u;
    for(U next=rest&a[u];next;next&=next-1)ends^=g[rest*n+__builtin_ctz(next)];
    g[m*n+u]=ends;
    endmask^=ends;
    if(__builtin_parity(ends))s[m]|=1u<<u;
    if(__builtin_popcount(m)%2){
     d[u]^=ends;
     for(U last=ends;last;last&=last-1)d[__builtin_ctz(last)]^=1u<<u;
    }
   }
   if(m!=full)boundary^=s[m];
   if(__builtin_popcount(m)%2==0)
    for(U left=full^m;left;left&=left-1)general[__builtin_ctz(left)]^=s[m]^endmask;
   if(n%2==0 && __builtin_popcount(m)%2==0 && __builtin_popcount(m)<=int(n)-4)
    for(U left=full^m;left;left&=left-1)h[__builtin_ctz(left)]^=s[m];
  }
  bool dz=true,hz=true;
  if(n%2==0 && h!=general)throw std::runtime_error("short support vs full endpoint sum");
  if(n>=3)for(U v=0;v<n;++v)
   if((__builtin_parity(d[v])^__builtin_parity(a[v]&general[v]))!=1)
    throw std::runtime_error("all-order path partition");
  for(U v=0;v<n;++v){dz &= d[v]==0;hz &= h[v]==0;}
  if(n%2){
   U total=0;
   for(U v=0;v<n;++v)total^=__builtin_parity(d[v]&~((1u<<(v+1))-1));
   if(total!=1)throw std::runtime_error("odd cover obstruction");
  }else if(n>=4){
   for(U v=0;v<n;++v)
    if((__builtin_parity(d[v])^__builtin_parity(a[v]&h[v]))!=1)
     throw std::runtime_error("even incompatibility");
  }
  auto &c=counts[n];++c[0];c[1]+=boundary==0;c[2]+=dz;
  if(n%2==0){c[3]+=hz;c[4]+=dz&&hz;c[5]+=boundary==0&&dz&&hz;}
 }
 for(U n=2;n<=9;++n){std::cout<<n;for(auto x:counts[n])std::cout<<" "<<x;std::cout<<"\n";}
}
