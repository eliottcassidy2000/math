// Exact anchored Wilson search for fixed upper bounds on the three center scores.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
using C=array<int,3>;
constexpr int N=53;
struct E{int y; C a; array<U,3> owner{};};
vector<E> B[N];
int reference[N],pick[N],usage[3][56],lev;
U initial,root;
uint64_t nodes,budget;
bool truncated;
U bit(int i){return 1ULL<<i;}
int pop(U x){return __builtin_popcountll(x);}
void change(const E& e,int d){for(int j=0;j<3;j++)usage[j][e.a[j]]+=d;}
struct Try{const E* e; U t,force; C q;};
bool dfs(U open,U active,U keep,C left){
 if(++nodes>budget){truncated=true;return false;}
 if(!open)return true;
 int ed=lev-pop(initial&~active), h=-1;size_t mn=999;
 vector<Try> L; C unavoidable{};
 for(U S=open;S;S&=S-1){
  int j=__builtin_ctzll(S);vector<Try> tmp; C lb{9,9,9};
  for(const auto& e:B[j]){
   Try p{&e,0,0,{}}; C imm{};bool ok=true;
   for(int c=0;c<3;c++){
    U u=e.owner[c]&active;
    int w=usage[c][e.a[c]];
    imm[c]=w+bool(u&keep);
    if(imm[c]>left[c]){ok=false;break;}
    if(u&&w==left[c])p.force|=u;
    p.q[c]=w+bool(u);p.t|=u;
   }
   if(ok&&pop(p.force)<=ed){
    tmp.push_back(p);
    for(int c=0;c<3;c++)lb[c]=min(lb[c],imm[c]);
   }
  }
  if(tmp.empty())return false;
  for(int c=0;c<3;c++)if((unavoidable[c]+=lb[c])>left[c])return false;
  if(tmp.size()<mn){mn=tmp.size();h=j;L=move(tmp);}
 }
 for(const auto& p:L){
  const E& e=*p.e;
  U optional=p.t&~(keep|p.force);
  for(U S=optional;;S=(S-1)&optional){
   U drop=S|p.force;
   if(pop(drop)<=ed){
    C next=left;
    // Changed owners are uncounted; all other collisions have distinct charges.
    for(int c=0;c<3;c++)next[c]-=p.q[c]-bool(e.owner[c]&drop);
    change(e,1);pick[h]=e.y;
    if(dfs((open^bit(h))|drop,active^drop,keep|(p.t&~drop),next))return true;
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
 if(n!=N||argc!=7)return 2;
 C cap{stoi(argv[1]),stoi(argv[2]),stoi(argv[3])};
 int start=stoi(argv[4]),stop=stoi(argv[5]);budget=stoull(argv[6]);
 if(*min_element(cap.begin(),cap.end())<0||start<0||stop>50||start>stop)return 2;
 int owner[3][56];for(auto&r:owner)fill(r,r+56,-1);
 E marker{};
 for(int h=0;h<n;h++){
  int k;cin>>reference[h]>>k;pick[h]=reference[h];
  while(k--){
   E e{};cin>>e.y;for(int& a:e.a)cin>>a;
   if(!h){marker=e;continue;}
   if(!equal(e.a.begin(),e.a.end(),marker.a.begin(),[](int p,int q){return p!=q;}))continue;
   if(e.y==reference[h]){initial|=bit(h);for(int c=0;c<3;c++){
     if(owner[c][e.a[c]]!=-1)return 3;
     owner[c][e.a[c]]=h;
   }}
   B[h].push_back(e);
  }
  if(reference[h]<0)root|=bit(h);
 }
 if(pop(root)!=2||pop(initial)!=50)return 3;
 for(int h=1;h<N;h++)for(auto&e:B[h])for(int c=0;c<3;c++)
  if(owner[c][e.a[c]]>=0)e.owner[c]=bit(owner[c][e.a[c]]);
 for(lev=start;lev<=stop;lev++){
  truncated=false;nodes=0;
  bool yes=dfs(root,initial,0,cap);
  cerr<<cap[0]<<cap[1]<<cap[2]<<" edits "<<lev<<" nodes "<<nodes<<" "
      <<(yes?"witness":truncated?"capped":"excluded")<<endl;
  if(yes){for(int y:pick)cout<<y<<" ";cout<<endl;return 0;}
  if(truncated)return 1;
 }
}
