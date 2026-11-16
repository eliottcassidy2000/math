#include <array>
#include <iostream>
#include <random>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <cstring>
using namespace std;
using Rows=array<int,12>;
constexpr int MAX=1<<12;
struct DP{
 int a[MAX][12]{};
 int run(Rows const&r,int n){
  int N=1<<n;
  memset(a,0,N*12*sizeof(int));
  for(int u=0;u<n;u++)a[1<<u][u]=1;
  for(int s=1;s<N-1;s++) for(int e=s;e;e&=e-1) {
    int i=__builtin_ctz(e),h=a[s][i];
    if(!h)continue;
    for(int e2=r[i]&~s;e2;e2&=e2-1){
      int j=__builtin_ctz(e2);
      a[s|(1<<j)][j]+=h;
    }
  }
  int h=0;for(int i=0;i<n;i++)h+=a[N-1][i];return h;
 }
 int paths(int i,int j,DP const&b,int n)const{
   int s0=1<<i,mask=((1<<n)-1)^s0^(1<<j),sum=0,full=(1<<n)-1;
   for(int s=mask;;s=(s-1)&mask){
     int t=s|s0;
     sum+=a[t][i]*b.a[full^t][j];
     if(s==0)break;
   }
   return sum;
 }
};
Rows conv(Rows r,int n){for(int i=0;i<n;i++)r[i]^=((1<<n)-1)^(1<<i);return r;}
void print(Rows const&r,int n){for(int i=0;i<n;i++)cout<<r[i]<<' ';cout<<'\n';}
void flip(Rows&r,int i,int j){r[i]^=1<<j;r[j]^=1<<i;}
int main(int argc,char**argv){
  static DP A,B;
  int n=argc>1?atoi(argv[1]):11,P=argc>2?atoi(argv[2]):30011,limit=argc>3?atoi(argv[3]):500000;
  if(n<2||n>12||limit<1||P<1||P>479001598)return 1;
  mt19937_64 gen(15163152523ull+n);
  vector<pair<int,int>> edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.emplace_back(i,j);
  int hits=0; bool seen=false;
  for(int it=1;it<=limit;it++){
    Rows R={}; auto k=gen(); int bits=64;
    for(auto [i,j]:edges){if(!bits){k=gen();bits=64;} if(k&1) R[i]|=1<<j;else R[j]|=1<<i;k>>=1;--bits;}
    int h=A.run(R,n);
    if(h!=P&&h!=P+2)continue;
    hits++;
    if(!seen){cout<<"count "<<h<<" iteration "<<it<<'\n';print(R,n);seen=true;}
    int t=2*P+2-h;
    for(int e1=-1;e1<(int)edges.size();e1++){
      Rows r=R;int h0=h;
      if(e1>=0){auto [i,j]=edges[e1];flip(r,i,j);h0=A.run(r,n);}
      else A.run(r,n);
      B.run(conv(r,n),n);
      for(int e2=e1+1;e2<(int)edges.size();e2++){
        auto [i,j]=edges[e2];if(r[j]>>i&1)swap(i,j);
        int d=A.paths(j,i,B,n)-A.paths(i,j,B,n);
        if(h0+d==t){
          flip(r,i,j);
          if(A.run(r,n)!=t)return 2;
          cout<<"FOUND "<<h<<" -> "<<t<<" iteration "<<it<<" hits "<<hits<<" e1 "<<e1<<" e2 "<<e2<<" edge(s) ";
          if(e1>=0)cout<<edges[e1].first<<','<<edges[e1].second<<' ';
          cout<<i<<','<<j<<'\n';print(R,n);print(r,n);return 0;
        }
      }
    }
  }
  cout<<"No pair found: hits "<<hits<<'\n';
}
