// Independent n8 ZERO-WEDGE delta audit: S7-orbit BFS via adjacent transpositions, then
// vertex-subset matching DP and wedge-first disjoint-edge counting.
// No producer headers, kernels, representatives or receipts are loaded.
#include <array>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <numeric>
#include <algorithm>
using namespace std;
void need(bool x,const char* s){if(!x)throw runtime_error(s);}
int pop(unsigned x){return __builtin_popcount(x);}
struct Counts{int m,m2,m3,m4,w,p21,residual;};
Counts evaluate(const array<unsigned,8>&a){
 array<array<int,5>,256> dp{};dp[0][0]=1;
 for(unsigned s=1;s<256;s++){
  unsigned vb=s&-s,t=s^vb;int v=__builtin_ctz(vb);dp[s]=dp[t];
  for(unsigned nb=a[v]&t;nb;nb&=nb-1){unsigned ub=nb&-nb;for(int j=1;j<=4;j++)dp[s][j]+=dp[t^ub][j-1];}
 }
 int w=0,p21=0;
 for(int v=0;v<8;v++){
  int d=pop(a[v]);w+=d*(d-1)/2;
  for(unsigned nn=a[v];nn;nn&=nn-1){unsigned ub=nn&-nn;
   for(unsigned kk=nn^ub;kk;kk&=kk-1){unsigned wb=kk&-kk;p21+=dp[255u^((1u<<v)|ub|wb)][1];}
  }
 }
 auto p=dp[255];need(w+p[2]==p[1]*(p[1]-1)/2,"edge-pair partition");
 int z=41*p[1]-32*p[2]+8*p[3]-4*p[4]+8*p21;
 return {p[1],p[2],p[3],p[4],w,p21,z};
}
int main(){
 vector<pair<int,int>> e7,e8;int idx7[7][7]{};int idx8[8][8]{};
 for(int i=0;i<7;i++)for(int j=i+1;j<7;j++){idx7[i][j]=idx7[j][i]=int(e7.size());e7.emplace_back(i,j);}
 for(int i=0;i<8;i++)for(int j=i+1;j<8;j++){idx8[i][j]=idx8[j][i]=int(e8.size());e8.emplace_back(i,j);}
 array<vector<pair<unsigned,unsigned>>,6> swaps;
 for(int v=0;v<6;v++)for(int x=0;x<7;x++)if(x!=v&&x!=v+1)swaps[v].emplace_back(1u<<idx7[v][x],1u<<idx7[v+1][x]);
 vector<uint8_t> seen(1u<<21);vector<pair<unsigned,unsigned>> representatives;uint64_t covered=0;
 for(unsigned g=0;g<(1u<<21);g++)if(!seen[g]){
  vector<unsigned> queue{g};seen[g]=1;
  for(size_t i=0;i<queue.size();i++)for(const auto& pairs:swaps){
   unsigned x=queue[i];for(auto [p,q]:pairs)if(bool(x&p)!=bool(x&q))x^=p|q;
   if(!seen[x]){seen[x]=1;queue.push_back(x);}
  }
  need(5040%queue.size()==0,"orbit does not divide S7 order");
  representatives.emplace_back(g,unsigned(queue.size()));covered+=queue.size();
 }
 need(covered==(1u<<21),"n7 coverage incomplete");need(representatives.size()==1044,"wrong n7 orbit count");
 array<unsigned,5041> orbit_hist{};uint64_t cases=0,zeros=0,digest=1469598103934665603ULL;int minimum=1000000,maximum=-1000000;
 vector<unsigned> zero_masks;vector<Counts> zero_counts;
 for(auto [g,os]:representatives){
  orbit_hist[os]++;
  array<unsigned,8>base{};unsigned base8=0;
  for(int k=0;k<21;k++)if(g>>k&1){auto[u,v]=e7[k];base[u]|=1u<<v;base[v]|=1u<<u;base8|=1u<<idx8[u][v];}
  for(unsigned attach=0;attach<128;attach++){
   auto a=base;unsigned encoded=base8;a[7]=attach;
   for(int i=0;i<7;i++)if(attach>>i&1){a[i]|=128u;encoded|=1u<<idx8[i][7];}
   Counts c=evaluate(a);minimum=min(minimum,c.residual);maximum=max(maximum,c.residual);
   need(c.residual>=0,"negative n8 graph residual");
   if(c.residual==0){zeros++;zero_masks.push_back(encoded);zero_counts.push_back(c);}
   digest^=uint64_t(c.residual);digest*=1099511628211ULL;cases++;
  }
 }
 need(cases==133632&&zeros==2&&minimum==0,"certificate summary mismatch");
 cout<<"{\"method\":\"adjacent-transposition orbit BFS; vertex-subset matching DP; wedge-first p21\",\"n\":8,\"seven_representatives\":"<<representatives.size()<<",\"seven_covered\":"<<covered<<",\"augmentations\":"<<cases<<",\"minimum\":"<<minimum<<",\"maximum\":"<<maximum<<",\"zero_augmentations\":"<<zeros<<",\"residual_digest\":"<<digest<<",\"orbit_size_histogram\":{";
 bool sep=false;for(int z=1;z<=5040;z++)if(orbit_hist[z]){if(sep)cout<<',';sep=true;cout<<'"'<<z<<"\":"<<orbit_hist[z];}
 cout<<"},\"zero_witnesses\":[";
 for(size_t i=0;i<zero_masks.size();i++){if(i)cout<<',';const auto&c=zero_counts[i];cout<<"{\"mask\":"<<zero_masks[i]<<",\"m\":"<<c.m<<",\"m2\":"<<c.m2<<",\"m3\":"<<c.m3<<",\"m4\":"<<c.m4<<",\"w\":"<<c.w<<",\"p21\":"<<c.p21<<'}';}
 cout<<"],\"status\":\"PASS\"}\n";
}
