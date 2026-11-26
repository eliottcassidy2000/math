// Find low collision assignments by branching on which reference owners
// are displaced by each new contact. Retained owners become permanent.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
constexpr int N=53;
struct E{int y;array<int,3> a;U owners;};
vector<E> B[N]; int reference[N],pick[N],q[3][56],lev,cap;
U initial,root;
uint64_t nodes,budget=30000000;
bool truncated;
U bit(int i){return 1ULL<<i;}
int pop(U x){return __builtin_popcountll(x);}
int used(const E& e){return q[0][e.a[0]]+q[1][e.a[1]]+q[2][e.a[2]];}
void change(const E& e,int d){for(int i=0;i<3;i++)q[i][e.a[i]]+=d;}
bool dfs(U open,U active,U keep,int left){
 if(++nodes>budget){truncated=true;return false;}
 if(!open)return true;
 int ed=lev-pop(initial&~active), h=-1;unsigned mn=999;
 vector<const E*> L;
 for(U S=open;S;S&=S-1){
  int j=__builtin_ctzll(S);vector<const E*> tmp;
  for(const auto& e:B[j]){
   U T=e.owners&active;
   if(used(e)+pop(T)-min(ed,pop(T&~keep))<=left)tmp.push_back(&e);
  }
  if(tmp.empty())return false;
  if(tmp.size()<mn){mn=tmp.size();h=j;L=move(tmp);}
 }
 for(auto p:L){
  const E&e=*p;U T=e.owners&active,M=T&~keep;
  int val=used(e)+pop(T);
  for(U S=M;;S=(S-1)&M){
   const int c=val-pop(S);
   if(c<=left&&pop(S)<=ed){
    change(e,1);pick[h]=e.y;
    if(dfs((open^bit(h))|S, active^S, keep|(T&~S),left-c))return true;
    change(e,-1);pick[h]=reference[h];
    if(truncated)return false;
   }
   if(!S)break;
  }
 }
 return false;
}
int main(int argc,char**argv){
 int n;cin>>n;
 if(n!=N||argc<4)return 2;
 cap=stoi(argv[1]);int start=stoi(argv[2]),stop=stoi(argv[3]);
 if(argc>4)budget=stoull(argv[4]);
 int owner[3][56];
 for(auto&r:owner)fill(r,r+56,-1);
 E marker{};
 for(int h=0;h<n;h++){
  int k;cin>>reference[h]>>k;pick[h]=reference[h];
  while(k--){
   E e{};cin>>e.y;for(int& a:e.a)cin>>a;
   if(!h){marker=e;continue;}
   if(!equal(e.a.begin(),e.a.end(),marker.a.begin(),[](int p,int q){return p!=q;}))continue;
   if(e.y==reference[h]){initial|=bit(h);for(int i=0;i<3;i++)owner[i][e.a[i]]=h;}
   B[h].push_back(e);
  }
  if(reference[h]<0)root|=bit(h);
 }
 if(pop(root)!=2||pop(initial)!=50)return 3;
 for(int h=1;h<N;h++){
  // A displaced reference is never useful at this same branch: an
  // already fixed edge hits it, so reverting would only add collisions.
  // Keep it available anyway to make the bounded search exact on ties.
  for(auto&e:B[h])for(int j=0;j<3;j++)if(owner[j][e.a[j]]>=0)e.owners|=bit(owner[j][e.a[j]]);
 }
 for(lev=start;lev<=stop;lev++){
  truncated=false;nodes=0;
  bool yes=dfs(root,initial,0,cap);
  cerr<<"edits "<<lev<<" cap "<<cap<<" nodes "<<nodes<<" result "<<(yes?"witness":truncated?"capped":"excluded")<<endl;
  if(yes){for(int y:pick)cout<<y<<" ";cout<<endl;return 0;}
  if(truncated)return 1;
 }
}
