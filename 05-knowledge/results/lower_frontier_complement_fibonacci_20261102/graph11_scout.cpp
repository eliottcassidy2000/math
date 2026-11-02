// Bounded probe only. This is NOT an exhaustive graph certificate.
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace std;
struct Result {int f,m,w,p21; array<int,6> matching;};
Result evaluate(uint64_t bits){
 int a[11]{},d[11]{},e=0;
 for(int i=0;i<11;i++)for(int j=i+1;j<11;j++,e++)if(bits>>e&1){a[i]|=1<<j;a[j]|=1<<i;}
 int W=0,p21=0;for(int i=0;i<11;i++){d[i]=__builtin_popcount((unsigned)a[i]);W+=d[i]*(d[i]-1)/2;}
 array<array<int,6>,2048>match{};match[0][0]=1;
 for(int s=1;s<2048;s++){
  int i=__builtin_ctz((unsigned)s),t=s^(1<<i);match[s]=match[t];
  for(int q=a[i]&t;q;q&=q-1){int j=__builtin_ctz((unsigned)q);for(int k=1;k<=5;k++)match[s][k]+=match[t^(1<<j)][k-1];}
 }
 for(int i=0;i<11;i++)for(int j=i+1;j<11;j++)if(a[i]>>j&1)
  for(int v=0;v<11;v++)if(v!=i&&v!=j){int z=d[v]-((a[v]>>i)&1)-((a[v]>>j)&1);p21+=z*(z-1)/2;}
 auto c=match[2047];int f=499*c[1]-320*c[2]+80*c[3]-20*c[4]+5*c[5]+80*p21+56*W;
 return {f,c[1],W,p21,c};
}
int main(){
 mt19937_64 rng(110740);int minimum=0;uint64_t witness=0;int sampled=0,negative=0;
 for(int t=0;t<100000;t++){
  uint64_t g=0;int density=1+t%19;
  for(int e=0;e<55;e++)if(rng()%20<density)g|=1ull<<e;
  auto r=evaluate(g);sampled++;if(r.f<0)negative++;
  if(r.f<minimum){minimum=r.f;witness=g;}
 }
 auto z=evaluate(witness);
 cout<<"{\"scope\":\"seeded probe, not exhaustive\",\"graphs\":"<<sampled
 <<",\"negative_samples\":"<<negative<<",\"minimum_including_empty\":"<<minimum
 <<",\"witness_bits\":"<<witness<<",\"matching_counts\":[";
 for(int k=0;k<=5;k++){if(k)cout<<",";cout<<z.matching[k];}
 cout<<"],\"wedges\":"<<z.w<<",\"p21\":"<<z.p21<<"}\n";
}
