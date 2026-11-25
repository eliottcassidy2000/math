#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using I=std::uint64_t;
using State=std::array<I,18>;
int ix(int a,int b,int c){return (3*a+b)*2+c;}
State push(State old,int w){
 State result{};
 for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<2;c++){
  I x=old[ix(a,b,c)];
  result[ix(w?2:1,b,c)]+=x;
  result[ix(a,b?b:(w?1:2),c)]+=x;
  result[ix(a,b,1)]+=x;
 }
 return result;
}
I best[2]={~I(0),~I(0)}, low[2]={0,0}, nearword[2]={0,0}, n=0, E=0;
int J[]={162,129}, S[]={11,16};
struct Reading{I height,L,K;};
Reading reading(State s,int chart){
 I Z=0,h=0,R=0,C=0;
 for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){
  I x=s[ix(a,b,c)];
  Z+=x;
  h+=(!c)*x;
  R+=(1+a)*x;
  C+=(!b)*x+(!c&&a)*x;
 }
 assert(R>=Z && C>=1);
 return {h*(111-S[chart])+29*(C-1)+(R-h)*J[chart]+S[chart],2*(320*(R-h)+h*(chart?191:158)),29*h};
}
I q(I x){
 I r=0;
 for(;x;x>>=1)if(x&1)r+=((x>>1)+1)>>1;
 return r;
}
I reverse(I x,int b){
 I y=0;
 for(int i=0;i<b;i++,x>>=1)y=(y<<1)|(x&1);
 return y;
}
void solve(Reading x,I T,I word,int c){
 if(x.height>T)return;
 T-=x.height;
 const I L=x.L,K=x.K,lo=T/L+(T%L!=0),hi=T/(L-K);
 if(lo>hi)return;
 int b=0;I B=1;
 while(lo>=2*B && __uint128_t(K)*(lo/(2*B)-1)>__uint128_t(L+K)*(2*B-1))B*=2,++b;
 const auto eval=[&](I u){return __uint128_t(L)*u-2*__uint128_t(K)*q(u);};
 for(I p=lo/B;p<=hi/B;p++){
  I l=0,r=B;
  while(l<r){
   I m=(l+r)/2;
   if(eval(p*B+reverse(m,b))>T)l=m+1;else r=m;
  }
  if(l<B){
   const I u=p*B+reverse(l,b);
   if(lo<=u&&u<=hi&&eval(u)==T)
    std::cout<<"hit "<<c<<" "<<T+x.height<<" "<<word<<" "<<u<<std::endl;
  }
 }
}
void visit(State s,int length,int max,I word){
 if(length>=max)return;
 auto z=push(s,0);
 I w=word*2;
 ++n;
 for(int c=0;c<2;c++){
  auto x=reading(z,c); I H=x.height,T=E-1;
  I D=H>=T?H-T:T-H;
  if(D<best[c]){best[c]=D;low[c]=H;nearword[c]=w;}
  for(I t=E-3;t<=E+1;t+=2)solve(x,t,w,c);
 }
 visit(z,length+1,max,w);
 visit(push(s,1),length+1,max,w+1);
}
int main(int argc,char**argv){
 if(argc==1){std::cout<<"Supply exponents, for example 7 13 17 19\n";return 0;}
 for(int i=1;i<argc;i++){
  int p=std::stoi(argv[i]);assert(p>=3&&p<=20);
  E=(I(1)<<(p-1))*((I(1)<<p)-1);
  int L=1;
  // Only search for exact endpoint targets (upper query E+1).
  while(129*((I(1)<<L)-1)<=E+1) ++L;
  State s{};s[0]=1;n=0;best[0]=best[1]=~I(0);
  visit(push(s,1),1,L,1);
  std::cout<<"p "<<p<<" max length "<<L<<" strong words "<<n<<" near (within scanned lengths) ";
  for(int c=0;c<2;c++)std::cout<<low[c]<<":"<<nearword[c]<<" ";
  std::cout<<"\n";
 }
}