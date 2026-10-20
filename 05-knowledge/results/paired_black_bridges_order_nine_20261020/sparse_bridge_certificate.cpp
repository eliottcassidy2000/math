#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
struct DSU{vector<int>p,s;DSU(int n):p(n),s(n,1){iota(p.begin(),p.end(),0);}int f(int x){return p[x]==x?x:p[x]=f(p[x]);}bool add(int x,int y){x=f(x);y=f(y);if(x==y)return false;if(s[x]<s[y])swap(x,y);p[y]=x;s[x]+=s[y];return true;}};
struct Edge{int a,b;uint32_t t;int forest;};
vector<vector<pair<int,int>>>g;vector<int>tin,lo;vector<Edge>es;vector<int>bridges;int timer=0;
void dfs(int x,int from){tin[x]=lo[x]=++timer;for(auto [y,e]:g[x]){if(e==from)continue;if(tin[y])lo[x]=min(lo[x],tin[y]);else{dfs(y,e);lo[x]=min(lo[x],lo[y]);if(lo[y]>tin[x])bridges.push_back(e);}}}
int main(int argc,char**argv){int n=argc>1?stoi(argv[1]):9;string prefix=argc>2?argv[2]:"census_n"+to_string(n);int m=(n-1)*(n-2)/2;uint32_t full=(1u<<m)-1;vector<int>map(1u<<m);ifstream f(prefix+".map",ios::binary);f.read((char*)map.data(),map.size()*sizeof(int));assert(f.gcount()==(streamsize)(map.size()*sizeof(int)));int N=*max_element(map.begin(),map.end())+1;DSU a(N),b(N);int ix[16][16];vector<pair<int,int>>co;for(int y=1;y<n-1;y++)for(int x=n;x>y+1;x--){ix[x][y]=co.size();co.push_back({x,y});}vector<int>sig;for(auto[x,y]:co)sig.push_back(ix[n+1-y][n+1-x]);int s=m/2;vector<uint32_t>l(1u<<s),h(1u<<(m-s));for(uint32_t t=1;t<l.size();t++){int k=__builtin_ctz(t);l[t]=l[t&(t-1)]|(1u<<sig[k]);}for(uint32_t t=1;t<h.size();t++){int k=__builtin_ctz(t);h[t]=h[t&(t-1)]|(1u<<sig[k+s]);}
 uint64_t total=0;for(uint32_t k=1u<<(m-1);k>0;k--){uint32_t t=k-1,rt=l[t&((1u<<s)-1)]|h[t>>s];if(t==rt||t>rt||t>(full^rt))continue;total++;int x=map[t],y=map[full^t];if(a.add(x,y))es.push_back({x,y,t,1});else if(b.add(x,y))es.push_back({x,y,t,2});}
 g.resize(N);for(int e=0;e<(int)es.size();e++){auto r=es[e];g[r.a].push_back({r.b,e});g[r.b].push_back({r.a,e});}tin.resize(N);lo.resize(N);int comp=0;for(int x=0;x<N;x++)if(!tin[x]){comp++;dfs(x,-1);}sort(bridges.begin(),bridges.end(),[](int i,int j){return es[i].t<es[j].t;});
 ofstream out(prefix+"_sparse_certificate.tsv");out<<"tiling\tu\tv\tforest\n";for(auto r:es)out<<r.t<<'\t'<<r.a<<'\t'<<r.b<<'\t'<<r.forest<<'\n';
 ofstream js(prefix+"_sparse_audit.json");js<<"{\"n\":"<<n<<",\"vertices\":"<<N<<",\"all_Q_edges_scanned\":"<<total<<",\"certificate_edges\":"<<es.size()<<",\"components\":"<<comp<<",\"bridges\":[";bool first=true;for(int e:bridges){if(!first)js<<',';first=false;auto r=es[e];js<<"{\"tiling\":"<<r.t<<",\"ends\":["<<r.a<<','<<r.b<<"]}";}js<<"]}\n";cerr<<"scanned="<<total<<" sparse="<<es.size()<<" bridges="<<bridges.size()<<endl;
}
