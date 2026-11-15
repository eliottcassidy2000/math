#include <bits/stdc++.h>
using namespace std;
using U=__uint128_t;
bool check(int M,int L){
 int D=1<<L, m=M/(M&-M);
 if(M*m>=128)throw 1;
 U goal=(U(1)<<(M*m))-1;
 vector<unsigned char> p(2*M*D);
 for(int x=1;x<int(p.size());++x){int u=x/2;p[x]=(p[u]+(x&1?0:(u+1)/2))%M;}
 vector<U> hit_for_x(M*m);
 for(int i=0;i<M*m;++i)hit_for_x[i]=U(1)<<i;
 for(int s=0;s<D;++s){
  int lo=s+1>=D-s-1?0:s+1, hi=lo==0?s+1:D, c=lo==0?s:s+D;
  for(int u=0;u<2*M*D;u+=D)for(int v=0;v<2*M*D;v+=D){
    U h0=0,h1=0;
    for(int b=lo;b<hi;++b){
      int k=(u+b)%m, f=(p[u+b]+p[v+c-b])%M;
      h0|=hit_for_x[f*m+k];
      h1|=hit_for_x[((f+(b&1))%M)*m+k];
      if(h0==goal && h1==goal)break;
    }
    if(h0!=goal || h1!=goal)return false;
  }
 }
 return true;
}
int main(){
 for(auto [M,L]:vector<pair<int,int>>{{7,12},{9,12},{10,11},{12,12},{16,10}}){
   if(!check(M,L)) throw runtime_error("unsupported modulus/length");
   cout<<"fixed-sum certificate "<<M<<" "<<L<<" PASS"<<endl;
 }
}
