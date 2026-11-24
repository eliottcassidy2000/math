#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <array>
using namespace std;
using U=uint64_t;
const int n=51;
const U ALL=(1ull<<n)-1;
long double L[n+1];
struct State {array<U,n> A; array<int,n> start,end;};
long double go(State &S,U rows,U cols,int depth) {
 U b=rows;int imin=-1,dmin=n;
 long double Lsum=1;
 while(b){
  int i=__builtin_ctzll(b);b&=b-1;
  int d=__builtin_popcountll(S.A[i]&cols);
  if(!d)return 0;
  if(!depth)Lsum*=L[d];
  else if(d<dmin)imin=i,dmin=d;
 }
 if(!depth)return Lsum;
 long double total=0;
 int i=imin,x=S.start[i];
 U a=S.A[i]&cols&~(1ull<<x);
 while(a){
  int j=__builtin_ctzll(a);a&=a-1;
  int y=S.end[j],u=S.start[y],v=S.end[x];
  U old=S.A[y];
  S.A[y]&=~(1ull<<x);S.start[y]=x;S.end[x]=y;
  total+=go(S,rows^(1ull<<i),cols^(1ull<<j),depth-1);
  S.A[y]=old;S.start[y]=u;S.end[x]=v;
 }
 return total;
}
int main(int argc,char**argv) {
 int depth=atoi(argv[1]); if(depth<0||depth>=n-1)return 1;
 for(int i=1;i<=n;i++){U F;cin>>F;L[i]=ldexpl((long double)F,-40);}
 State S{};
 for(int i=0;i<n;i++){cin>>S.A[i];S.start[i]=S.end[i]=i;}
 // All operations have positive normal results; this dyadic inflation safely
 // exceeds the accumulated IEEE error with up to 51 additions per depth.
 long double upper=go(S,ALL,ALL,depth)*(1+ldexpl(1.L,-30));
 printf("%.12Lf %La\n",log2l(upper),upper);
}
