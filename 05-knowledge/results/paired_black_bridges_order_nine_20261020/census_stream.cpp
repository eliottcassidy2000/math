#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

struct OnlineBridges {
 vector<int> e,cc,sz,par,last; vector<uint32_t> pe; int it=0;
 set<uint32_t> bridges;
 OnlineBridges(int n):e(n),cc(n),sz(n,1),par(n,-1),last(n,0),pe(n,0){iota(e.begin(),e.end(),0);iota(cc.begin(),cc.end(),0);}
 int f2(int v){if(v==-1)return -1; return e[v]==v?v:e[v]=f2(e[v]);}
 int fc(int v){v=f2(v); return cc[v]==v?v:cc[v]=fc(cc[v]);}
 void root(int v){v=f2(v);int r=v,child=-1;uint32_t old=0;while(v!=-1){int p=f2(par[v]);uint32_t edge=pe[v];par[v]=child;pe[v]=old;cc[v]=r;child=v;v=p;old=edge;}sz[r]=sz[child];}
 void merge(int a,int b){++it;vector<int>x,y;int l=-1;while(l==-1){if(a!=-1){a=f2(a);x.push_back(a);if(last[a]==it){l=a;break;}last[a]=it;a=par[a];}if(b!=-1){b=f2(b);y.push_back(b);if(last[b]==it){l=b;break;}last[b]=it;b=par[b];}}for(int v:x){e[v]=l;if(v==l)break;assert(bridges.erase(pe[v])==1);}for(int v:y){e[v]=l;if(v==l)break;assert(bridges.erase(pe[v])==1);}}
 void add(int a,int b,uint32_t id){a=f2(a);b=f2(b);if(a==b)return;int ca=fc(a),cb=fc(b);if(ca!=cb){if(sz[ca]>sz[cb]){swap(a,b);swap(ca,cb);}root(a);par[a]=b;pe[a]=id;cc[a]=b;sz[cb]+=sz[a];bridges.insert(id);}else merge(a,b);}
};

struct Node {uint32_t rep;uint64_t mass=0,blue=0,H=0;bool sc=false,strong=false;int c3=0;vector<int> scores;};
int n,m;uint32_t full;vector<int> cmap;vector<uint32_t> rlow,rhigh;int split;int adj[16],pathv[16],bp[16][16];vector<Node> nodes;
uint32_t refl(uint32_t t){return rlow[t&((1u<<split)-1)]|rhigh[t>>split];}
void dfs(int used,int depth,uint32_t mask,int u){
 if(depth==n){Node &r=nodes[u];r.H++;if(mask==refl(r.rep))r.sc=true;if(cmap[mask]<0){uint32_t rm=refl(mask);cmap[mask]=u;r.mass++;if(rm==mask)r.blue++;else if(cmap[rm]<0){cmap[rm]=u;r.mass++;}else assert(cmap[rm]==u);}else assert(cmap[mask]==u);return;}
 int poss=((depth?adj[pathv[depth-1]]:(1<<n)-1)&~used);
 while(poss){int v=__builtin_ctz((unsigned)poss);poss&=poss-1;uint32_t next=mask;for(int i=0;i<depth-1;++i)if(adj[v]&(1<<pathv[i]))next|=1u<<bp[i][depth];pathv[depth]=v;dfs(used|(1<<v),depth+1,next,u);}
}
int main(int argc,char**argv){
 n=argc>1?stoi(argv[1]):9;string prefix=argc>2?argv[2]:"census_n"+to_string(n);assert(n>=3&&n<=9);auto begin=chrono::steady_clock::now();
 m=(n-1)*(n-2)/2;full=(1u<<m)-1;vector<pair<int,int>> coords;int ix[16][16];for(auto &row:ix)fill(std::begin(row),std::end(row),-1);for(int y=1;y<n-1;y++)for(int x=n;x>y+1;x--){ix[x][y]=coords.size();coords.push_back({x,y});bp[n-x][n-y]=coords.size()-1;}
 vector<int>sig(m);for(int i=0;i<m;i++){auto [x,y]=coords[i];sig[i]=ix[n+1-y][n+1-x];}
 split=m/2;rlow.resize(1u<<split);rhigh.resize(1u<<(m-split));for(uint32_t t=1;t<rlow.size();t++){int k=__builtin_ctz(t);rlow[t]=rlow[t&(t-1)]|(1u<<sig[k]);}for(uint32_t t=1;t<rhigh.size();t++){int k=__builtin_ctz(t);rhigh[t]=rhigh[t&(t-1)]|(1u<<sig[k+split]);}
 cmap.assign(1u<<m,-1);uint64_t assigned=0;
 for(uint32_t t=0;t<=full;t++){
  if(cmap[t]>=0)continue;
  int u=nodes.size();Node node;node.rep=t;
  fill(adj,adj+n,0);for(int i=0;i<n-1;i++)adj[i]|=1<<(i+1);for(int k=0;k<m;k++){auto [x,y]=coords[k];int i=n-x,j=n-y;if(t&(1u<<k))adj[j]|=1<<i;else adj[i]|=1<<j;}
  int reach[16];for(int i=0;i<n;i++){reach[i]=adj[i]|(1<<i);int d=__builtin_popcount((unsigned)adj[i]);node.scores.push_back(d);node.c3-=d*(d-1)/2;}node.c3+=n*(n-1)*(n-2)/6;sort(node.scores.begin(),node.scores.end());
  for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]&(1<<k))reach[i]|=reach[k];
  node.strong=true;for(int i=0;i<n;i++)if(reach[i]!=(1<<n)-1)node.strong=false;
  nodes.push_back(node);dfs(0,0,0,u);auto&r=nodes[u];assert(r.H%2==1);assert(r.sc?(r.mass%2==1&&r.blue%2==1):(r.mass%4==2&&r.blue==0));uint64_t f=r.sc?r.mass:r.mass/2;assert(r.H%f==0&&(r.H/f)%2==1);assigned+=r.mass;
  if(u%1000==0)cerr<<"classes "<<u+1<<" assigned "<<assigned<<" last "<<t<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<endl;
 }
 assert(assigned==1u<<m);cerr<<"classification done nodes "<<nodes.size()<<endl;
 ofstream mapout(prefix+".map",ios::binary);mapout.write((char*)cmap.data(),cmap.size()*sizeof(int));mapout.close();
 OnlineBridges g(nodes.size()),bgraph(nodes.size());vector<uint64_t> deg(nodes.size(),0),bdeg(nodes.size(),0);uint64_t edges=0,loops=0,be=0,bl=0;
 for(uint32_t t=0;t<(1u<<(m-1));t++){
  uint32_t rt=refl(t);int a=cmap[t],b=cmap[full^t];
  if(rt==t){be++;if(a==b)bl++;bdeg[a]++;bdeg[b]++;bgraph.add(a,b,t);continue;}
  if(t>rt||t>(full^rt))continue;
  ++edges;if(a==b)loops++;deg[a]++;deg[b]++;g.add(a,b,t);
 }
 int comps=0,bcomps=0,active=0,sc=0,classes=0;for(int u=0;u<(int)nodes.size();u++){assert(deg[u]==(nodes[u].mass-nodes[u].blue)/2);assert(bdeg[u]==nodes[u].blue);if(g.fc(u)==u)comps++;if(bgraph.fc(u)==u)bcomps++;if(deg[u])active++;sc+=nodes[u].sc;classes+=nodes[u].sc?1:2;}
 ofstream out(prefix+".json");out<<"{\n\"n\":"<<n<<",\"tilings\":"<<(1u<<m)<<",\"classes\":"<<classes<<",\"merged\":"<<nodes.size()<<",\"SC\":"<<sc<<",\"paired_black\":{\"edges\":"<<edges<<",\"loops\":"<<loops<<",\"components\":"<<comps<<",\"active_vertices\":"<<active<<",\"cycle_rank\":"<<edges-nodes.size()+comps<<",\"bridges\":[";
 bool first=true;for(uint32_t t:g.bridges){if(!first)out<<",";first=false;out<<"{\"tiling\":"<<t<<",\"reflection\":"<<refl(t)<<",\"complement\":"<<(full^t)<<",\"ends\":["<<cmap[t]<<","<<cmap[full^t]<<"]}";}out<<"]},\"blue\":{\"edges\":"<<be<<",\"loops\":"<<bl<<",\"components\":"<<bcomps<<",\"bridge_count\":"<<bgraph.bridges.size()<<"},\n\"nodes\":[\n";
 first=true;for(int u=0;u<(int)nodes.size();u++){auto&r=nodes[u];if(!first)out<<",\n";first=false;out<<"{\"node\":"<<u<<",\"representative\":"<<r.rep<<",\"mass\":"<<r.mass<<",\"blue\":"<<r.blue<<",\"H\":"<<r.H<<",\"aut\":"<<r.H/(r.sc?r.mass:r.mass/2)<<",\"SC\":"<<(r.sc?"true":"false")<<",\"strong\":"<<(r.strong?"true":"false")<<",\"degree\":"<<deg[u]<<",\"c3\":"<<r.c3<<",\"scores\":[";for(int i=0;i<n;i++){if(i)out<<",";out<<r.scores[i];}out<<"]}";}
 out<<"\n],\"elapsed_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<"}\n";cerr<<"DONE n="<<n<<" nodes="<<nodes.size()<<" bridges="<<g.bridges.size()<<" elapsed="<<chrono::duration<double>(chrono::steady_clock::now()-begin).count()<<endl;
}
