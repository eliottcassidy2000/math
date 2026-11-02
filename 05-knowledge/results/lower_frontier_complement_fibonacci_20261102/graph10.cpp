// Exact integer graph certificate. No external graph census is used.
// The seven-vertex orbit construction extends Session138's graph9 producer.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;
using U=uint64_t;
struct Graph {
 int n; vector<pair<int,int>> edges; int ix[10][10]{};
 explicit Graph(int n):n(n){for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
  ix[i][j]=ix[j][i]=edges.size();edges.push_back({i,j});}}
 vector<int> adjacency(U g) const {vector<int>a(n);for(int e=0;e<(int)edges.size();e++)if(g>>e&1){auto[i,j]=edges[e];a[i]|=1<<j;a[j]|=1<<i;}return a;}
 U augment(U g,int neighbors,const Graph& larger) const {
  U h=0;for(int e=0;e<(int)edges.size();e++)if(g>>e&1){auto[i,j]=edges[e];h|=1ull<<larger.ix[i][j];}
  for(int i=0;i<n;i++)if(neighbors>>i&1)h|=1ull<<larger.ix[i][n];return h;
 }
 U canon(U g) const {
  auto a=adjacency(g);vector<int>color(n);for(int i=0;i<n;i++)color[i]=__builtin_popcount((unsigned)a[i]);
  // Canonical equitable refinement: previous colour is included in the signature.
  for(int round=0;round<n;round++){
   vector<vector<int>>sig(n);for(int i=0;i<n;i++){
    sig[i].push_back(color[i]);vector<int>neigh;for(int j=0;j<n;j++)if(a[i]>>j&1)neigh.push_back(color[j]);
    sort(neigh.begin(),neigh.end());sig[i].insert(sig[i].end(),neigh.begin(),neigh.end());}
   auto labels=sig;sort(labels.begin(),labels.end());labels.erase(unique(labels.begin(),labels.end()),labels.end());
   vector<int>next(n);for(int i=0;i<n;i++)next[i]=lower_bound(labels.begin(),labels.end(),sig[i])-labels.begin();
   bool same=true;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if((next[i]==next[j])!=(color[i]==color[j]))same=false;
   color=next;if(same)break;
  }
  vector<int>p(n);iota(p.begin(),p.end(),0);sort(p.begin(),p.end(),[&](int i,int j){return color[i]!=color[j]?color[i]<color[j]:i<j;});
  vector<pair<int,int>>blocks;for(int i=0;i<n;){int j=i+1;while(j<n&&color[p[i]]==color[p[j]])j++;blocks.push_back({i,j});i=j;}
  U best=~U(0);auto rec=[&](auto&&self,int b)->void{
   if(b==(int)blocks.size()){U h=0;for(int e=0;e<(int)edges.size();e++){auto[i,j]=edges[e];h|=((g>>ix[p[i]][p[j]])&1)<<e;}best=min(best,h);return;}
   auto[l,r]=blocks[b];do{self(self,b+1);}while(next_permutation(p.begin()+l,p.begin()+r));};rec(rec,0);return best;
 }
};
struct Data{vector<array<int,5>> match;vector<int>a,d;int m,W,p21,F;};
Data data(const Graph&G,U g){
 Data z;int n=G.n;z.a=G.adjacency(g);z.d.resize(n);z.m=__builtin_popcountll(g);z.W=0;z.p21=0;
 for(int i=0;i<n;i++){z.d[i]=__builtin_popcount((unsigned)z.a[i]);z.W+=z.d[i]*(z.d[i]-1)/2;}
 z.match.resize(1<<n);z.match[0][0]=1;
 for(int s=1;s<(1<<n);s++){
  int i=__builtin_ctz((unsigned)s),t=s^(1<<i);z.match[s]=z.match[t];
  for(int ns=z.a[i]&t;ns;ns&=ns-1){int j=__builtin_ctz((unsigned)ns);for(int k=1;k<=4;k++)z.match[s][k]+=z.match[t^(1<<j)][k-1];}
 }
 for(auto[i,j]:G.edges)if(z.a[i]>>j&1)for(int v=0;v<n;v++)if(v!=i&&v!=j){
  int d=z.d[v]-((z.a[v]>>i)&1)-((z.a[v]>>j)&1);z.p21+=d*(d-1)/2;}
 auto c=z.match.back();z.F=25*c[1]-16*c[2]+4*c[3]-c[4]+4*z.p21;return z;
}
vector<int> extensions(const Graph&G,U g){
 auto z=data(G,g);int n=G.n,full=(1<<n)-1;vector<int>L(n),out(1<<n);int Q[10][10]{};
 for(int i=0;i<n;i++){
  int nd=0;for(int j=0;j<n;j++)if(z.a[i]>>j&1)nd+=z.d[j];
  int Wdel=z.W-z.d[i]*(z.d[i]-1)/2-nd+z.d[i];
  int ell=Wdel+z.d[i]*(z.m-z.d[i]+1)-nd;
  auto c=z.match[full^(1<<i)];L[i]=25-16*c[1]+4*c[2]-c[3]+4*ell;
  for(int j=0;j<n;j++)if(i!=j){Q[i][j]=z.m-z.d[i]-z.d[j]+((z.a[i]>>j)&1);if(Q[i][j]<0)throw runtime_error("negative outside edge count");}
 }
 out[0]=z.F;for(int s=1;s<=full;s++){int i=__builtin_ctz((unsigned)s),t=s^(1<<i);int v=out[t]+L[i];for(int ns=t;ns;ns&=ns-1)v+=4*Q[i][__builtin_ctz((unsigned)ns)];out[s]=v;}return out;
}
int main(int argc,char**argv){try{
 Graph g7(7),g8(8),g9(9),g10(10);mt19937_64 rng(8102710);
 int direct_checks=0;
 for(int t=0;t<200;t++){
  U g=rng()&((1ull<<36)-1);auto v=extensions(g9,g);
  for(int s=0;s<512;s+=17){if(v[s]!=data(g10,g9.augment(g,s,g10)).F)throw runtime_error("extension identity");direct_checks++;}
  vector<int>p(9);iota(p.begin(),p.end(),0);shuffle(p.begin(),p.end(),rng);U h=0;
  for(int e=0;e<36;e++)if(g>>e&1){auto[i,j]=g9.edges[e];h|=1ull<<g9.ix[p[i]][p[j]];}
  if(g9.canon(g)!=g9.canon(h))throw runtime_error("canonical relabelling");
 }
 if(argc>1&&string(argv[1])=="--scout"){
  int minimum=0;U witness=0;for(int t=0;t<20000;t++){
   U g=0;unsigned threshold=(t%19)+1;for(int e=0;e<45;e++)if(rng()%20<threshold)g|=1ull<<e;
   int f=data(g10,g).F;if(f<minimum){minimum=f;witness=g;}}
  cout<<"{\"scout_only\":true,\"graphs\":20000,\"direct_extension_checks\":"<<direct_checks<<",\"minimum\":"<<minimum<<",\"witness\":"<<witness<<"}\n";return 0;
 }
 vector<array<uint32_t,21>>perms;array<int,7>p={0,1,2,3,4,5,6};do{array<uint32_t,21>a;for(int e=0;e<21;e++){auto[i,j]=g7.edges[e];a[e]=1u<<g7.ix[p[i]][p[j]];}perms.push_back(a);}while(next_permutation(p.begin(),p.end()));
 vector<uint8_t>seen(1<<21);vector<U>r7;uint64_t covered=0;
 for(uint32_t g=0;g<(1u<<21);g++)if(!seen[g]){
  r7.push_back(g);for(auto&a:perms){uint32_t h=0;for(uint32_t b=g;b;b&=b-1)h|=a[__builtin_ctz(b)];if(!seen[h]){seen[h]=1;covered++;}}}
 if(covered!=(1u<<21))throw runtime_error("coverage7");
 unordered_set<U>r8set;for(U g:r7)for(int s=0;s<128;s++)r8set.insert(g8.canon(g7.augment(g,s,g8)));
 vector<U>r8(r8set.begin(),r8set.end());sort(r8.begin(),r8.end());cerr<<"covered7="<<covered<<" reps8="<<r8.size()<<"\n";
 unordered_set<U>r9set;uint64_t aug9=0;
 for(U g:r8){for(int s=0;s<256;s++){r9set.insert(g9.canon(g8.augment(g,s,g9)));aug9++;}if(aug9%262144==0)cerr<<"aug9="<<aug9<<" reps9="<<r9set.size()<<"\n";}
 vector<U>r9(r9set.begin(),r9set.end());sort(r9.begin(),r9.end());cerr<<"reps9="<<r9.size()<<"\n";
 uint64_t cases=0,zeros=0;int minimum=100000000;U witness=0;map<pair<int,int>,uint64_t>zero_stats;
 for(U g:r9){auto vals=extensions(g9,g);for(int s=0;s<512;s++){
  int f=vals[s];cases++;if(f<minimum){minimum=f;witness=g9.augment(g,s,g10);}if(f==0){zeros++;U h=g9.augment(g,s,g10);auto z=data(g10,h);zero_stats[{z.m,z.W}]++;}
  if(f<0){cerr<<"negative witness="<<witness<<" value="<<minimum<<"\n";}}
  if(cases%16777216==0)cerr<<"cases10="<<cases<<" min="<<minimum<<"\n";
 }
 cout<<"{\"n\":10,\"covered7\":"<<covered<<",\"reps7\":"<<r7.size()<<",\"reps8\":"<<r8.size()<<",\"augmentations9\":"<<aug9<<",\"reps9\":"<<r9.size()<<",\"augmentations10\":"<<cases<<",\"minimum\":"<<minimum<<",\"zeros\":"<<zeros<<",\"witness\":"<<witness<<",\"direct_extension_checks\":"<<direct_checks<<",\"zero_stats\":[";
 bool first=true;for(auto [key,count]:zero_stats){if(!first)cout<<",";first=false;cout<<"{\"edges\":"<<key.first<<",\"wedges\":"<<key.second<<",\"augmentations\":"<<count<<"}";}cout<<"]}\n";
 if(minimum<0)return 2;
}catch(exception&e){cerr<<e.what()<<"\n";return 1;}}
