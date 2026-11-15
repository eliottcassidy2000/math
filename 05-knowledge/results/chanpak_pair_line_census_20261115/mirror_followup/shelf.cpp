#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;
uint32_t rev(uint32_t x){
 x=((x&0x55555555U)<<1)|((x>>1)&0x55555555U);
 x=((x&0x33333333U)<<2)|((x>>2)&0x33333333U);
 x=((x&0x0f0f0f0fU)<<4)|((x>>4)&0x0f0f0f0fU);
 x=((x&0x00ff00ffU)<<8)|((x>>8)&0x00ff00ffU);
 return (x<<16)|(x>>16);
}
int main(int argc,char**argv){
 int lim=argc==2?atoi(argv[1]):16;if(lim<3||lim>25)return 2;
 vector<uint32_t> I(1U<<lim);
 for(uint32_t x=1;x<I.size();++x){
    uint32_t q=x>>1;
    I[x]=I[q]+((x&1)==0?((q+1)>>1):0);
 }
 for(int n=3;n<=lim;++n){
    const uint32_t N=1U<<n;
    vector<uint64_t> bucket(N);
    for(uint32_t x=0;x<N;++x){
      auto y=rev(x)>>(32-n);
      int e=int(x&1)-int(x>>(n-1));
      uint64_t id=3ULL*(N+x-y)+(e+1);
      bucket[x]=(id<<n)|I[x];
    }
    sort(bucket.begin(),bucket.end());
    unsigned r=1,best=1;
    for(size_t i=1;i<bucket.size();++i){
      uint64_t a=bucket[i-1],b=bucket[i];
      if(a==b)continue;
      if((a>>n)==(b>>n) && b==a+1)best=max(best,++r);
      else r=1;
    }
    cout<<n<<" "<<best<<endl;
 }
}
