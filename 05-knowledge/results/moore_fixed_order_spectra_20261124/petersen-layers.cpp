#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
constexpr int n=10, B=n/2, N=1<<n;
using Poly=array<int,n>;
bool graph[n][n];
int bound[N];
uint64_t mask(const array<int,n>& q) {
 uint64_t k=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
  k=2*k+graph[q[i]][q[j]];
 return k;
}
array<int,n> row(const array<int,n>& q) {
 array<int,n> a{};
 for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
  if(graph[q[i]][q[j]]) a[j]|=1<<i; else a[i]|=1<<j;
 return a;
}
Poly paths(const array<int,n>& out, int split) {
 Poly dp[N][n]{};
 for(int i=0;i<n;i++) dp[1<<i][i][0]=1;
 for(int a=1;a<N-1;a++) for(int ends=a;ends;ends&=ends-1) {
  int i=__builtin_ctz(ends);
  for(int nxt=out[i]&~a;nxt;nxt&=nxt-1) {
   int j=__builtin_ctz(nxt),t=((split>>i)&1)==((split>>j)&1);
   for(int z=0;z<=bound[a]-1;z++) dp[a|(1<<j)][j][z+t]+=dp[a][i][z];
  }
 }
 Poly h{};for(int i=0;i<n;i++)for(int z=0;z<n;z++) h[z]+=dp[N-1][i][z];
 return h;
}
int sum(const Poly&h) {int s=0;for(int x:h)s+=x;return s;}
void show(const Poly&h) {for(int j=0;j<n;j++)cout<<(j?",":"")<<h[j];}
int main() {
 for(int i=0;i<5;i++)for(auto [j,k]: {pair{i,(i+1)%5},pair{i,5+i},pair{5+i,5+(i+2)%5}})
  graph[j][k]=graph[k][j]=true;
 for(int i=0;i<N;i++)bound[i]=__builtin_popcount((unsigned)i);
 array<int,n> q{}; for(int i=0;i<n;i++)q[i]=i;
 set<uint64_t> seen;
 map<int,int> freq,alt;
 map<array<int,4>,array<int,3>> cuttypes;
 int raw=0, examples=0;
 do {
  bool ok=true;for(int i=0;i<n-1;i++)if(graph[q[i]][q[i+1]]){ok=false;break;}
  if(!ok)continue;raw++;
  auto k=mask(q);if(!seen.insert(k).second)continue;
  auto a=row(q);
  Poly h=paths(a,(1<<B)-1);
  int cross=0, matches=0, direct[2]{};
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)if(graph[q[i]][q[j]]) {
    if(i<B&&j>=B)cross++;else direct[i/B]++;
  }
  array<int,B> p{};for(int j=0;j<B;j++)p[j]=B+j;
  do {bool yes=true;for(int j=0;j<B;j++)yes&=graph[q[j]][q[p[j]]];matches+=yes;}while(next_permutation(p.begin(),p.end()));
  array<int,4> c{h[0],cross,matches,min(direct[0],direct[1])};
  cuttypes[c][0]++;
  for(int i=0;i<n-1;i++) {
   if(i==B-1)continue;
   swap(q[i],q[i+1]);
   bool good=(i==0||!graph[q[i-1]][q[i]])&&(i==n-2||!graph[q[i+1]][q[i+2]]);
   uint64_t l=mask(q);swap(q[i],q[i+1]);
   if(!good||l<k)continue;
   auto aa=a;aa[i]^=1<<(i+1);aa[i+1]^=1<<i;
   Poly g=paths(aa,(1<<B)-1);
   assert(g[0]==h[0]);
   int gap=abs(sum(g)-sum(h));freq[gap]++; cuttypes[c][1]++;cuttypes[c][2]+=(gap==2);
   if(gap==2) {
    alt[h[0]]++;
    if(examples++<12) {
      cout<<"seam "<<k<<" "<<l<<" @"<<i<<" heights "<<sum(h)<<" "<<sum(g)<<" spectra ";
      show(h);cout<<" / ";show(g);cout<<" order ";
      for(int v:q)cout<<v;cout<<"\n";
    }
   }
  }
 }while(next_permutation(q.begin(),q.end()));
 int total=0;for(auto [d,f]:freq)total+=f;
 assert(raw==95520 && seen.size()==796 && total==1330 && freq[0]==18 && freq[2]==10);
 assert((alt==map<int,int>{{156,10}}));
 cout<<"raw "<<raw<<" masks "<<seen.size()<<" moves "<<total<<" zero "<<freq[0]<<" two "<<freq[2]<<"\n";
 for(auto [c,v]:cuttypes){
   cout<<"cut";for(int x:c)cout<<" "<<x;cout<<" : ";for(int x:v)cout<<x<<" ";cout<<"\n";
 }
 cout<<"seams by alternating count";for(auto [a,c]:alt)cout<<" "<<a<<":"<<c;cout<<"\n";
}