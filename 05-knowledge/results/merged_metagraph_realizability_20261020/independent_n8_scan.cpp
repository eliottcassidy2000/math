// Independent direct-arc reimplementation of the complete n=8 scan.
// It has no precomputed move masks, affine bit transforms, atlas, or isomorphism code.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
struct UnionFind {
 vector<uint32_t>p;uint32_t count;
 UnionFind(uint32_t N):p(N),count(N){iota(p.begin(),p.end(),0);}
 uint32_t find(uint32_t v){while(v!=p[v]){p[v]=p[p[v]];v=p[v];}return v;}
 void merge(uint32_t a,uint32_t b){a=find(a);b=find(b);if(a!=b){if(a>b)swap(a,b);p[b]=a;--count;}}
};
int main(){
 const int n=8,m=21;const uint32_t N=1u<<m;vector<pair<int,int>>f;int index[n][n];for(auto&r:index)for(auto&v:r)v=-1;
 for(int i=0;i<n;++i)for(int j=i+2;j<n;++j){index[i][j]=f.size();f.push_back({i,j});}
 vector<array<int,n>>words;
 for(int a=0;a<n;++a)for(int b=a+2;b<n;++b){vector<int>w(n);iota(w.begin(),w.end(),0);int t=w[a];w.erase(w.begin()+a);w.insert(w.begin()+b,t);array<int,n>q;copy(w.begin(),w.end(),q.begin());words.push_back(q);}
 UnionFind u(N);uint64_t moves=0;
 for(uint32_t x=0;x<N;++x){
  bool A[n][n]={};for(int i=0;i<n-1;++i)A[i][i+1]=true;
  for(int k=0;k<m;++k){auto[i,j]=f[k];A[i][j]=(x>>k)&1u;A[j][i]=!A[i][j];}
  for(const auto&p:words){
   bool ok=true;for(int i=0;i<n-1;++i)if(!A[p[i]][p[i+1]]){ok=false;break;}if(!ok)continue;
   ++moves;uint32_t y=0;for(int k=0;k<m;++k){auto[i,j]=f[k];if(A[p[i]][p[j]])y|=1u<<k;}u.merge(x,y);
  }
 }
 uint32_t ordinary=u.count;
 for(uint32_t x=0;x<N;++x){uint32_t y=0;for(int k=0;k<m;++k){auto[i,j]=f[k];if((x>>k)&1u)y|=1u<<index[n-1-j][n-1-i];}u.merge(x,y);}
 cout<<"{\"n\":8,\"tilings\":"<<N<<",\"valid_moves\":"<<moves<<",\"ordinary_components\":"<<ordinary<<",\"merged_components\":"<<u.count<<"}"<<endl;
 if(moves!=8912896||ordinary!=6880||u.count!=3528)return 1;
}
