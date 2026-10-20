#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
// Floating scout only. Every reported candidate requires an exact replay.
int main(int argc,char**argv){
 uint64_t lim=argc>1?std::strtoull(argv[1],nullptr,10):100000001;
 if(lim>(uint64_t(1)<<34))throw std::runtime_error("declared source ceiling");
 long double best=-1,maxD=0; uint64_t steps=0,maxdepth=0,maxSource=0,maxTarget=0,firstTwo=0;
 std::cout<<std::setprecision(22);
 for(uint64_t n=3;n<=lim;n+=2){
  uint64_t x=n,k=0,m=0;long double D=0;
  do {
   if(x>(UINT64_MAX-1)/3)throw std::runtime_error("overflow");
   uint64_t v=3*x+1;int a=__builtin_ctzll(v);x=v>>a;k+=a;++m;++steps;
   D=std::ldexp(3*D,-a)+(a>=3?1-std::ldexp(4.0L,-a):0);
   if(m>10000)throw std::runtime_error("depth cap");
  }while(x>=n);
  long double ratio=D/(n-x);
  if(ratio>best){best=ratio;std::cout<<"record "<<n<<" "<<x<<" "<<m<<" "<<k<<" "<<D<<" "<<ratio<<"\n";}
  if(D>maxD){maxD=D;maxSource=n;maxTarget=x;}
  if(D>=2 && firstTwo==0){firstTwo=n;std::cout<<"first_D_at_least_two "<<n<<" "<<x<<" "<<D<<"\n";}
  if(m>maxdepth)maxdepth=m;
  if(ratio>=1){std::cout<<"candidate "<<n<<"\n";return 2;}
 }
 std::cout<<"limit "<<lim<<" steps "<<steps<<" maxdepth "<<maxdepth<<" maxD "<<maxD<<" source "<<maxSource<<" target "<<maxTarget<<"\n";
}
