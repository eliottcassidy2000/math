#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
using namespace std;
constexpr int N=3250;
using Mask=bitset<N>;
int n,d,s;vector<int>E[N];Mask A[N],R[N],C[N];
void r2(int u){R[u]=A[u];R[u].set(u);for(int x:E[u])R[u]|=A[x];}
void r3(int u){C[u]=R[u];for(int x:E[u])C[u]|=R[x];}
void add(int u,int v){E[u].push_back(v);E[v].push_back(u);A[u][v]=A[v][u]=1;}
int addr(int i,int j){return d+1+s*i+j;}
bool pack(int b,int c,mt19937 &g,vector<pair<int,int>>&p,bool rnd){
 vector<int>L,T;
 for(int i=0;i<s;i++){int u=addr(b,i);if(none_of(E[u].begin(),E[u].end(),[&](int w){return w>d&&(w-d-1)/s==c;}))L.push_back(u);}
 for(int j=0;j<s;j++){int v=addr(c,j);if(none_of(E[v].begin(),E[v].end(),[&](int w){return w>d&&(w-d-1)/s==b;}))T.push_back(v);}
 assert(L.size()==T.size()); if(L.empty())return false;
 vector<vector<int>> h(L.size());
 for(int i=0;i<(int)L.size();i++)for(int j=0;j<(int)T.size();j++)if(!C[L[i]][T[j]])h[i].push_back(j);
 if(rnd)for(auto &a:h)shuffle(a.begin(),a.end(),g);
 vector<int> z(L.size(),-1); vector<char> seen(L.size());
 auto aug=[&](auto &&self,int i)->bool{
   for(int j:h[i])if(!seen[j]){seen[j]=1;if(z[j]<0||self(self,z[j])){z[j]=i;return true;}}
   return false;
 };
 for(int i=0;i<(int)L.size();i++){fill(seen.begin(),seen.end(),0);if(!aug(aug,i))return false;}
 p.clear();for(int j=0;j<(int)T.size();j++)p.emplace_back(L[z[j]],T[j]);
 return true;
}
int main(int argc,char**argv){
 assert(argc==3);
 cin>>n;
 for(int u=0;u<n;u++){int m;cin>>m;while(m--){int v;cin>>v;if(v>u)add(u,v);}}
 d=E[0].size();s=d-1;assert(n<=N&&n==d*d+1);
 mt19937 g(stoi(argv[1]));int lim=stoi(argv[2]);
 for(int u=0;u<n;u++)r2(u);
 for(int u=0;u<n;u++)r3(u);
 vector<pair<int,int>> cand,p;
 for(int b=0;b<d;b++)for(int c=b+1;c<d;c++)cand.emplace_back(b,c);
 shuffle(cand.begin(),cand.end(),g);int done=0,edge=0;
 for(auto [b,c]:cand){
  if(!pack(b,c,g,p,true))continue;
  Mask aff;
  for(auto [u,v]:p){assert(!C[u][v]);aff|=R[u]|R[v];add(u,v);}
  for(auto [u,v]:p)aff|=A[u]|A[v];
  vector<int> Q;for(int i=0;i<n;i++)if(aff[i]){Q.push_back(i);r2(i);}
  for(int i:Q)r3(i);
  // Any newly appearing path of length <=3 has all endpoints in aff.
  done++;edge+=p.size();
  cerr<<"board "<<b<<" "<<c<<" edges "<<p.size()<<" totals "<<done<<" "<<edge<<"\n";
  if(done==lim)break;
 }
 cout<<n<<"\n";
 for(auto &x:E){cout<<x.size();sort(x.begin(),x.end());for(int v:x)cout<<" "<<v;cout<<"\n";}
}
