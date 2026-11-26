// Bounded, exact-at-each-node search on the residue/phase factorization.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>
using namespace std;
struct Choice {array<int,3> r,eta;int mask,j;};
vector<Choice> bank[52];
int initial[84],maxcap[84],pick[52];
uint64_t visits,limit;
bool phase;
struct DSU {
 uint8_t p[53],z[53],w[53];
 DSU(){for(int i=0;i<53;i++){p[i]=i;z[i]=0;w[i]=1;}}
 int root(int x,int& a)const {a=0;while(p[x]!=x){a^=z[x];x=p[x];}return x;}
 bool link(int x,int y,int c){
  int a,b; x=root(x,a);y=root(y,b);c^=a^b;
  if(x==y)return c==0;
  if(w[x]>w[y])swap(x,y);
  p[x]=y;z[x]=c;w[y]+=w[x];return true;
 }
};
struct State {
 array<int8_t,84> cap,first,sgn;
 DSU d;
 uint64_t live=(1ULL<<52)-1;
};
bool accept(State &S,int i,const Choice &e){
 if(phase&&e.mask!=3&&!S.d.link(i,52,e.mask/2))return false;
 for(int k=0;k<3;k++){
  int h=e.r[k],f=e.eta[k];
  if(!S.cap[h])return false;
  if(phase && ((maxcap[h]==1 && !S.d.link(i,52,f^(initial[h]/2)))||
               (maxcap[h]==2 && S.cap[h]==1 && !S.d.link(i,S.first[h],1^f^S.sgn[h]))))return false;
  --S.cap[h];S.first[h]=i;S.sgn[h]=f;
 }
 S.live^=1ULL<<i;
 return true;
}
bool dfs(const State &S){
 if(++visits>limit)return false;
 if(!S.live)return true;
 pair<int,int> top[1064];int L=0;
 int scores[136]{};
 for(uint64_t u=S.live;u;u&=u-1){
  int i=__builtin_ctzll(u),sc=0;
  for(int j=0;j<(int)bank[i].size();j++){
   const auto&e=bank[i][j];
   auto N=S;
   if(!accept(N,i,e))continue;
   top[L++]={i,j};++sc;
   for(int h:e.r)++scores[52+h];
  }
  if(!sc)return false;
  scores[i]=sc;
 }
 int m=5000,h=-1;
 for(uint64_t u=S.live;u;u&=u-1){int i=__builtin_ctzll(u);
  if(scores[i]<m){m=scores[i];h=i;}
 }
 for(int i=0;i<84;i++)if(S.cap[i]){
  int sc=scores[52+i];
  if(sc<S.cap[i])return false;
  if((sc+S.cap[i]-1)/S.cap[i]<m){m=(sc+S.cap[i]-1)/S.cap[i];h=52+i;}
 }
 for(int j=0;j<L && visits<=limit;j++){
  auto [i,l]=top[j]; const auto&e=bank[i][l];
  if(h<52?i!=h:find(e.r.begin(),e.r.end(),h-52)==e.r.end())continue;
  auto N=S; if(!accept(N,i,e))return false;
  pick[i]=e.j;
  if(dfs(N))return true;
 }
 return false;
}
int main(int argc,char**argv){
 phase=argc>=2&&string(argv[1])=="phase";
 limit=argc>2?stoull(argv[2]):100000;
 mt19937 gen(argc>3?stoi(argv[3]):14285);
 State S{};
 for(int i=0;i<84;i++){cin>>initial[i];S.cap[i]=maxcap[i]=__builtin_popcount(initial[i]);}
 for(auto &B:bank){
  int k,j=0;cin>>k;
  while(k--){
   Choice e;int y;cin>>y;
   for(int g=0;g<3;g++){int s;cin>>s;e.r[g]=28*g+s%28;e.eta[g]=s/28;}
   cin>>e.mask;e.j=j++;B.push_back(e);
  }
  shuffle(B.begin(),B.end(),gen);
 }
 if(!dfs(S)){cerr<<"No result in "<<visits<<" visited nodes.\n";return 1;}
 cerr<<"Profile after "<<visits<<" nodes.\n";
 for(int p:pick)cout<<p<<" ";
 cout<<endl;
}
