#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
#include <bit>
#include <cmath>
using namespace std;
using u64=uint64_t; using i64=int64_t;using u128=__uint128_t;
constexpr unsigned S=1u<<16;
u64 Vcoef,Ucoef,Xcoef,hcoef,c0;
array<u64,S> qtab,revtab;
u64 qfast(u64 x){
  u64 head=x>>48, q=qtab[head];
  for(int j=32;j>=0;j-=16){
    u64 a=(x>>j)&65535, R=revtab[a];
    q+=qtab[a]+(head>>1)*R+(head&1)*((R+(R&1))>>1);
    head=(head<<16)|a;
  }
  return q;
}
u64 rev(u64 v,int b){
  u64 r=revtab[v&65535];
  if(b>16)r=(r<<16)|revtab[(v>>16)&65535];
  return r>>(b>16?32-b:16-b);
}
u64 isqrt(u64 x){u64 q=sqrtl((long double)x);while((q+1)<=x/(q+1))q++;while(q&&q>x/q)q--;return q;}
u64 ceildiv(u64 a,u64 b){return a/b+(a%b!=0);}
inline u64 eval(u64 L,u64 t){return L*t-2*hcoef*qfast(t);}
int main(int argc,char**argv){
 for(unsigned i=1;i<S;i++){
  qtab[i]=qtab[i>>1]+((i&1)?((i>>1)+1)/2:0);
  revtab[i]=(revtab[i>>1]>>1)|((i&1)<<15);
 }
 for(u64 x=1234567,i=0;i<20000;i++){
  x=x*11400714819323198485ULL+11;u64 a=x>>3,q=0;
  for(u64 v=a;v;v>>=1)if(v&1)q+=((v>>1)+1)/2;
  assert(qfast(a)==q);
 }
 const bool eight=argc>5&&string(argv[5])=="eight";
 Vcoef=eight?13536:382;Ucoef=eight?5382:316;Xcoef=eight?35232:1280;hcoef=eight?651:29;c0=eight?2313:111;
 const int p=argc>1?stoi(argv[1]):31;
 const int delta=argc>2?stoi(argv[2]):1;
 assert(4<=p&&p<=31&&-3<=delta&&delta<=1);
 const u64 H=(u64(1)<<(p-1))*((u64(1)<<p)-1)+delta;
 if(H<c0){cout<<"below minimum"<<endl;return 0;}
 const u64 top=isqrt((H-c0)/Xcoef);
 const u64 start=argc>3?stoull(argv[3]):0, stop=argc>4?min<u64>(top,stoull(argv[4])):top;
 u64 checked=0, prefixes=0,hits=0;
 for(u64 s=start;s<=stop;s++){
  const u64 qS=qfast(s);
  for(int j=0;j<2;j++){
   const u64 L=Xcoef*s+(j?Ucoef:Vcoef),D=c0+(j?Vcoef:Ucoef)*s-2*hcoef*qS;
   if(D>H)continue;
   const u64 T=H-D;
   const u64 lo=max(s,ceildiv(T,L)),hi=T/(L-hcoef);
   if(hi<lo)continue;
   checked++;
   const int residue=(T%hcoef);const bool unit=L%hcoef;
   if(!unit&&residue)continue;
   if(hi-lo<6){
    for(u64 t=lo;t<=hi;t++)if((t%hcoef)*(L%hcoef)%hcoef==u64(residue) && eval(L,t)==T){
      ++hits;cout<<"hit "<<s<<" "<<t<<" "<<j<<" "<<H<<endl;
    }
    continue;
   }
   // Use the candidate interval's actual lower end, not its dyadic width.
   int b=0;u64 B=1;
   while(lo>=2*B&&u128(hcoef)*(lo/(2*B)-1)>u128(L+hcoef)*(2*B-1)){B<<=1;++b;if(b>=31)abort();}
   for(u64 p=lo/B;p<=hi/B;p++){
    ++prefixes;u64 a=0,r=B;
    while(a<r){const u64 m=(a+r)/2,t=p*B+rev(m,b);
      if(eval(L,t)>T)a=m+1;else r=m;
    }
    if(a<B){
      const u64 t=p*B+rev(a,b);
      if(lo<=t&&t<=hi&&(t%hcoef)*(L%hcoef)%hcoef==u64(residue)&&eval(L,t)==T){
       ++hits;cout<<"hit "<<s<<" "<<t<<" "<<j<<" "<<H<<endl;
      }
    }
   }
  }
 }
 cout<<(eight?"eight ":"six ")<<"range "<<start<<" "<<stop<<" H "<<H<<" feasible "<<checked<<" prefixes "<<prefixes<<" hits "<<hits<<endl;
}