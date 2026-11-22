#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;
using V=uint64_t;
using I=__int128_t;
using U=__uint128_t;
using W=array<V,4>;
W wide(U x){return {V(x),V(x>>64),0,0};}
W times(W a,V b){
 U carry=0;
 for(int i=0;i<4;i++){U w=U(a[i])*b+carry;a[i]=w;carry=w>>64;}
 return a;
}
W sum(W a,W b){
 U carry=0;
 for(int i=0;i<4;i++){U w=U(a[i])+b[i]+carry;a[i]=w;carry=w>>64;}
 return a;
}
W mul(W a,U b){
 W c=times(a,V(b>>64));
 return sum(times(a,V(b)),{0,c[0],c[1],c[2]});
}
bool ge(W const&a,W const&b){
 for(int i=3;i>=0;i--)if(a[i]!=b[i])return a[i]>b[i];
 return true;
}
struct Chart{V h,q,c,t,e,r,k;};
const Chart C4={2465,2293760,1344270,712944,1740,986,423385};
const Chart C8={37265,144834560,84839910,57602544,24360,14906,33838015};
const I N61=(I(1)<<60)*((I(1)<<61)-1)-1;
I integer(string s){I v=0;for(char c:s)v=10*v+c-'0';return v;}
string out(I x){if(x<0)return "-"+out(-x);string z;do{z+='0'+x%10;x/=10;}while(x);reverse(z.begin(),z.end());return z;}
int mod(I x,int p){int r=x%p;return r<0?r+p:r;}
I floor(I a,I b){I q0=a/b;return q0-(a%b<0);}
I ceiling(I a,I b){return -floor(-a,b);}
struct Byte {
 array<V,256> F{},R{};
 Byte(){for(int i=1;i<256;++i){F[i]=F[i>>1]+(i&1?0:(i>>1)+((i>>1)&1));R[i]=(R[i>>1]>>1)|((i&1)<<7);}}
}const B;
V f(V x){V ans=0;while(x>>8){int l=x&255;x>>=8;ans+=B.F[l]+(255-B.R[l])*x+(l<128&&(x&1));}return ans+B.F[x];}
I ht(Chart c,V x,V y){return I(c.q)*x*y+I(c.t)*x+I(c.c)*y+c.r*(x&1)+c.e*(y&1)+c.k+I(c.h)*(f(x)+f(y));}
V inv(V a,V p){
 if(p==1)return 0;
 int64_t s=1,t=0;V b=p;
 while(b){V q=a/b,r=a%b;a=b;b=r;int64_t u=s-int64_t(q)*t;s=t;t=u;}
 return (s%int64_t(p)+p)%p;
}
struct Root{int n,j,k;}; // n=-1 means all classes
struct Mod {
 int p;vector<int> iv,sq;
 Mod(int p0):p(p0),iv(p),sq(p,-1){
  for(int j=0;j<p;j++){iv[j]=inv(j,p);sq[j*j%p]=j;} // entry iv[0] is never read
 }
 Root roots(int a,int b,int c)const{
  if(!a){if(b)return {1,mod(-c*iv[b],p),0};return {c?0:-1,0,0};}
  int d=sq[mod(b*b-4*a*c,p)];if(d<0)return {0,0,0};
  int t=iv[(2*a)%p];return {d?2:1,mod((-b+d)*t,p),mod((-b-d)*t,p)};
 }
};
pair<V,V> combine(pair<V,V> v,int a,Mod const& p){
 auto [r,m]=v;int g=m%p.p;
 int d=mod(I(a)-r,p.p);
 if(!g)return d?make_pair(V(0),V(0)):v;
 return {r+m*(d*p.iv[g]%p.p),m*p.p};
}
pair<V,V> slope(long double z,V den){
 V a0=0,b0=1,a1=1,b1=0;
 const long double target=z;
 for(int i=0;i<40;i++){
  V x=z;
  if(b1 && x>(den-b0)/b1){
   V l=(den-b0)/b1,A=a0+l*a1,B=b0+l*b1;
   if(fabsl((long double)A/B-target)<fabsl((long double)a1/b1-target))return {A,B};
   return {a1,b1};
  }
  V A=a0+x*a1,B=b0+x*b1;
  a0=a1;b0=b1;a1=A;b1=B;
  z-=x;if(z<1e-18L||B==den)break;
  z=1/z;
 }
 return {a1,b1};
}
struct Branch {
 Chart c; I target,R; bool sw;
 V o,p,s,u,Q,T,C,K;
 Mod pa,pb;
 V slabs=0,lines=0,reads=0,hits=0;
 Branch(Chart src,I n,int i,int j,bool swap):c(src),target(n),sw(swap),
 o(swap?j:i),p(swap?i:j),s(swap?2:10),u(swap?10:2),
 Q(20*c.q),T(s*((swap?c.c:c.t)+c.q*p)),C(u*((swap?c.t:c.c)+c.q*o)),
 K(c.k+c.t*i+c.c*j+c.r*(i&1)+c.e*j+c.q*i*j),pa(c.h==37265?257:17),pb(29){
  R=target-K;
 }
 void block(V X,V width){
  slabs++;
  I den0=I(Q)*X+C;
  V y0=(R-I(T)*X)/den0;
  if(I(y0)<I(X))return;
  long double der=((long double)Q*R+(long double)T*C)/((long double)den0*den0);
  auto [a,b]=slope(der,width);
  auto bound=[&](V d,I R1,V T1,V C1,bool lower){
   I pos=X+d,den=Q*pos+C1;
   I F=I(a)*d*den+b*(R1-I(T1)*pos-I(y0)*den);
   return lower?ceiling(F,den):floor(F,den);
  };
  I upper=max(bound(0,R,T,C,false),bound(width-1,R,T,C,false));
  I R1=R-I(c.h)*(o+p);
  V T1=T+c.h*s,C1=C+c.h*u;
  V dmin=0;
  if(a && width>1){
   const W cache=times(sum(times(wide(R1),Q),wide(U(T1)*C1)),b);
   long double der1=sqrtl((long double)b/a*((long double)Q*R1+(long double)T1*C1));
   long double may=(der1-C1)/Q-X-0.5L;
   dmin=clamp(may,0.L,(long double)width-1);
   auto nonneg=[&](V d){I D=I(Q)*(X+d)+C1;return ge(mul(times(wide(D),a),D+Q),cache);};
   while(dmin>0&&nonneg(dmin-1))--dmin;
   while(dmin+1<width&&!nonneg(dmin))++dmin;
  }else if(!a)dmin=width-1;
  I lower=bound(dmin,R1,T1,C1,true);
  if(lower>upper)return;
  V bi=inv(a,b);
  for(I z=lower;z<=upper;z++){
   ++lines;
   I J=I(b)*y0+z;
   I QA=-I(Q)*a,LB=I(b)*T+Q*J-I(a)*den0;
   auto constant=[&](int p){return mod(mod(b,p)*mod(T,p)*mod(X,p)+mod(den0,p)*mod(J,p)-mod(b,p)*mod(R,p),p);};
   auto U=pa.roots(mod(QA,pa.p),mod(LB,pa.p),constant(pa.p));
   auto Vv=pb.roots(mod(QA,pb.p),mod(LB,pb.p),constant(pb.p));
   pair<V,V> base={V((z%I(b)+b)%b)*bi%b,b};
   for(int i=0;i<abs(U.n);i++){
    auto v=U.n<0?base:combine(base,i?U.k:U.j,pa);if(!v.second)continue;
    for(int j=0;j<abs(Vv.n);j++){
     auto w=Vv.n<0?v:combine(v,j?Vv.k:Vv.j,pb);if(!w.second)continue;
     for(V d=w.first;d<width;d+=w.second){
      I Y=I(y0)+(z-I(a)*d)/b;
      if(Y<I(X+d+sw))continue;
      V x=s*(X+d)+o,y=u*Y+p;
      if(sw)std::swap(x,y);
      ++reads;
      if(ht(c,x,y)==target){++hits;cout<<"HIT "<<x<<" "<<y<<endl;}
     }
    }
   }
  }
 }
 void run(V low,V high,V width){
  auto exact=[&](V x){return I(x)*x<=R/Q;};
  V l=low,r=high;
  if(!exact(l))return;
  while(l<r){V m=l+(r-l+1)/2;if(exact(m))l=m;else r=m-1;}
  for(V X=low;X<=l;X+=width)block(X,min(width,l-X+1));
  cout<<"branch "<<o<<","<<p<<","<<sw<<" rows "<<(l-low+1)
      <<" slabs "<<slabs<<" lines "<<lines<<" reads "<<reads<<" hits "<<hits<<endl;
 }
};
int main(int argc,char**argv){
 if(argc<5)return 2;
 int n=atoi(argv[1]);Chart c=n==4?C4:C8;
 V low=integer(argv[2]),high=integer(argv[3]),width=integer(argv[4]);
 I N=argc>5?integer(argv[5]):N61;
 // 128-bit exact anchored line numerators with the following fixed scope.
 if((n!=4&&n!=8)||low<1000000000000ULL||high>300000000000000ULL||low>high||
    !width||width>65536||N>N61)return 2;
 for(int i=0;i<10;i++){
  if((N-c.k-I(c.t)*i-c.r*(i&1))%5)continue;
  for(int j=0;j<2;j++)for(bool sw:{false,true}){
   Branch b(c,N,i,j,sw);
   if(b.R>0)b.run(low,high,width);
  }
 }
}
