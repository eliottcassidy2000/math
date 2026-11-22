#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using I=__int128;
using U=unsigned __int128;
using namespace std;
I P(int k){return I(1)<<k;}
string out(U i){if(!i)return "0";string s;while(i){s+=char('0'+i%10);i/=10;}reverse(s.begin(),s.end());return s;}
int main(int argc,char**argv){
 const int limit=argc>1?stoi(argv[1]):5;
 const int absorder=argc>2?stoi(argv[2]):0;
 const I target=P(121)-P(60)-1;
 int n;long long c,h,k0,k1;string mask;
 uint64_t tests=0,suffixes=0;int bank=0;
 while(cin>>n>>c>>h>>k0>>k1>>mask){
  ++bank;
  const I ymax=(target-min(k0,k1))/c;
  int L=0;while(P(L)<=ymax)++L;
  const int need=absorder?absorder-7-n:0;
  if(need>L||need<0)continue;
  if(need)L=need;
  auto accept=[&](I A,I A0,I B0,I B1,I z,int b,int cc,int d,int e){
   // L(A)=A0 A + z; if A starts with a nonempty all-one run, A+1 is a power of two.
   if(z>=target || A0>target-z) return;
   if(((target-z)|A0)%2==0) {
    // The leading run gives odd A. Unequal 2-adic valuations of the sides reject.
    auto vn=[](I t){uint64_t t0=(uint64_t)t;return t0?__builtin_ctzll(t0):64+__builtin_ctzll((uint64_t)(t>>64));};
    if(vn(target-z)!=vn(A0))return;
   }
   ++tests;
   I t=target-z;
   if(t%A0)return;
   A=t/A0;
   if(((A+1)&A)!=0)return;
   I y=B0+B1*A;
   if(y>ymax)return;
   int wid=0;for(I k=y;k;k>>=1)++wid;
   if(need&&wid!=need)return;
   I f=0,x=y;while(x){I p=x>>1;if(!(x&1))f+=p+(p&1);x=p;}
   assert((y&1?k1:k0)+c*y+h*f==target);
   cout<<"HIT "<<mask<<" "<<out(y)<<" run suffix "<<b<<","<<cc<<","<<d<<","<<e<<"\n";
  };
  // Cases 1,2,3 were already exhausted in the source packet. New cases are 4,5.
  for(int b=1;b<L-2;++b)for(int j=1;b+j<L-1;++j)
   for(int d=1;b+j+d<L;++d)for(int e=0;b+j+d+e<L&&e<=L;++e){
    if(e&&limit<5)break;
    // x is the suffix-independent address after 1^a 0^b 1^j:
    // x = B1p*A+B0p, p(x)=(2^b-1)A+1.
    I B1p=P(b+j),B0p=P(j)-1;
    I B1=B1p<< (d+e),B0=(B0p<<(d+e))+P(e)-1;
    I F1=P(b)-1+(P(d)-1)*B1p,F0=2+(P(d)-1)*B0p;
    ++suffixes;
    accept(0,c*B1+h*F1,B0,B1,(e?k1:k0)+c*B0+h*F0,b,j,d,e);
   }
 }
 cerr<<"bank="<<bank<<", suffixes="<<suffixes<<", exact division tests="<<tests<<"\n";
}
