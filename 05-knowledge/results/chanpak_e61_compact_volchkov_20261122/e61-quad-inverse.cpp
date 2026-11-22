#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using I=__int128;using U=unsigned __int128;
using namespace std;
I P(int k){return I(1)<<k;}
const I N=P(121)-P(60)-1;
int n,W,R;long long c,h,k0,k1;string key; int min0,min1,lim1;
uint64_t tuples=0,hit=0;
void emit(int a,int b){
 ++hit;cout<<"HIT "<<key<<" first/last "<<a<<" "<<b<<"\n";
}
void solve(int M,I A,I B,I C,I D){
 ++tuples;int lo=1,hi=min(lim1,M-min1);if(lo>hi)return;
 int odd=R%2;I kw=odd?c:(I)c+h;
 I z=kw*A*P(M)+(odd?k1-c+h*D:k0+h*(D-B+1))-N;
 I u=h*(odd?C:C-A), v=odd?c*(B+1):kw*B;
 auto f=[&](int a)->I{return u*P(a)+v*P(M-a)+z;};
 auto find=[&](int l,int r,int s){
   if(s*f(l)>0||s*f(r)<0)return;
   while(l<r){int m=(l+r)/2;if(s*f(m)<0)l=m+1;else r=m;}
   if(f(l)==0)emit(l,M-l);
 };
 if(odd){find(lo,hi,1);return;}
 // u,v < 0, so forward differences strictly decrease.  Mid is the first
 // possible maximum; include the second equal point on the right side.
 int l=lo,r=hi;
 while(l<r){
   int m=(l+r)/2;
   if(u*P(m)-v*P(M-m-1)>0)l=m+1;
   else r=m;
 }
 find(lo,l,1);
 if(l<hi)find(l+1,hi,-1);
}
void walk(int i,int left,I A,I B,I C,I D){
 if(i==R-1){solve(left,A,B,C,D);return;}
 const int bound=i==1?min0:1;
 const int upper=left-(R-i-2)-min1-1;
 for(int j=bound;j<=upper;++j){
   I k=P(j);
   if(i%2)walk(i+1,left-j,k*A,k*B,C+(k-1)*A,D+(k-1)*B+1);
   else walk(i+1,left-j,k*A,k*(B+1)-1,C,D);
 }
}
int main(int argc,char**argv){
 int tot=argc>1?stoi(argv[1]):123,cap=argc>2?stoi(argv[2]):7,fl=argc>3?stoi(argv[3]):6;
 int seen=0;
 while(cin>>n>>c>>h>>k0>>k1>>key){
   W=tot-7-n;
   if(W<=0||c*P(W-1)+min(k0,k1)>N)continue;
   if((c+h)*(P(W)-1)+max(k0,k1)<N)continue;
   ++seen;
   I top=(N-min(k0,k1))/c, Y=P(W-1),room=N-min(k0,k1)-c*Y;
   min0=1;while(min0+1<W&&top<Y+(Y>>(min0+1)))++min0;
   int t=0;while(t+1<W&&(Y>>(t+1))>room/h)++t;
   lim1=top<Y+(Y>>1)?1:W;
   for(R=fl;R<=cap;++R){
     if(t&&R%2==0)continue;
     min1=max(1,t);
     // With leading length forced to one, no search in it is needed.
     walk(1,W,1,-1,0,0);
   }
 }
 cerr<<"order="<<tot<<" run window="<<fl<<".."<<cap<<" charts="<<seen<<" middle tuples="<<tuples<<" hits="<<hit<<"\n";
}
