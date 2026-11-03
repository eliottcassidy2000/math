#include <cstdint>
#include <iostream>
#include <vector>
int main() {
  constexpr int n=21;
  uint32_t out[n]{};
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   bool rev=(i==0&&j==7)||(i==5&&j==9)||(i==9&&j==14)||(i==14&&j==20);
   if(rev)out[j]|=1u<<i;else out[i]|=1u<<j;
  }
  const uint32_t M=1u<<n;
  std::vector<uint64_t> d((size_t)M*n);
  for(int i=0;i<n;i++) d[((size_t)1<<i)*n+i]=1;
  for(uint32_t s=1;s<M;s++){
   uint32_t x=s;
   while(x){int i=__builtin_ctz(x);x&=x-1;auto v=d[(size_t)s*n+i];if(!v)continue;
    uint32_t y=out[i]&~s&(M-1);while(y){int j=__builtin_ctz(y);y&=y-1;d[(size_t)(s|(1u<<j))*n+j]+=v;}
   }
  }
  uint64_t h=0;for(int i=0;i<n;i++)h+=d[(size_t)(M-1)*n+i];
  std::cout<<"21-vertex consecutive attachment H="<<h<<" expected="<<265*17*33<<"\n";
  return h==265*17*33?0:1;
}
