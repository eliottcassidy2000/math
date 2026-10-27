#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
int n; vector<pair<int,int>> es,chosen,witness; array<int,16> parent,sz,d; uint64_t counts[16], hist[40]; int maxL;
int root(int x){while(parent[x]!=x)x=parent[x];return x;}
void visit(int k){
 if(k==(int)es.size()){
  ++counts[chosen.size()]; auto s=d;sort(s.begin(),s.begin()+n);int l1=0;for(int i=0;i<n;++i)l1+=abs(s[i]-i);
  if(l1%2 || l1/2>n-2){cerr<<"VIOLATION";exit(2);} int L=l1/2; ++hist[L];
  if(L>maxL){maxL=L;witness=chosen;}
  return;
 }
 visit(k+1);auto [a,b]=es[k];int ra=root(a),rb=root(b); if(ra==rb)return;
 if(sz[ra]<sz[rb]) { swap(ra,rb); }
 parent[rb]=ra;sz[ra]+=sz[rb]; ++d[a];--d[b];chosen.emplace_back(a,b);
 visit(k+1);
 chosen.pop_back(); --d[a];++d[b];sz[ra]-=sz[rb];parent[rb]=rb;
}
int main(int argc,char**argv){int upto=argc>1?atoi(argv[1]):8;for(n=2;n<=upto;++n){es.clear();chosen.clear();witness.clear();fill(begin(counts),end(counts),0);fill(begin(hist),end(hist),0);maxL=-1;for(int i=0;i<n;++i){parent[i]=i;sz[i]=1;d[i]=i;for(int j=i+1;j<n;++j)es.emplace_back(i,j);}visit(0);cout<<"n="<<n<<" total="<<accumulate(begin(counts),end(counts),uint64_t(0))<<" maxL="<<maxL<<" counts_by_edges=";for(int i=0;i<n;++i)cout<<counts[i]<<',';cout<<" histogram=";for(int i=0;i<=maxL;++i)cout<<i<<':'<<hist[i]<<',';cout<<" witness=";for(auto[a,b]:witness)cout<<a<<'-'<<b<<',';cout<<'\n'<<flush;}}
