// All nearest projected repairs by anchored exact enumeration.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
constexpr U ALL=(1ULL<<52)-1;
U bit(int i){return 1ULL<<i;}
struct E{int y;array<int,3> r;};
vector<E> R[52];
array<int,52> reference,sol;
set<array<int,52>> ans;
U init,nodes,block[52][56];
int drop,lim;
bool fits(const E&e,const array<U,3>&used){
 for(int z=0;z<3;z++)if(z!=drop && (used[z]&bit(e.r[z])))return false;
 return true;
}
void dfs(U open,U active,array<U,3>used){
 ++nodes;
 if(!open){ans.insert(sol);return;}
 int h=-1, mn=100000, rem=lim-__builtin_popcountll(init^active);
 vector<const E*> opts;
 for(U f=open;f;f&=f-1){
  int i=__builtin_ctzll(f);
  vector<const E*> v;
  for(const E&e:R[i])if(fits(e,used)&&__builtin_popcountll(block[i][e.y]&active)<=rem)v.push_back(&e);
  if(v.empty())return;
  if((int)v.size()<mn){h=i;mn=v.size();opts=move(v);}
 }
 for(auto ep:opts){
  auto&e=*ep;U hit=block[h][e.y]&active;
  auto V=used;for(int z=0;z<3;z++)if(z!=drop)V[z]|=bit(e.r[z]);
  sol[h]=e.y;
  dfs((open^bit(h))|hit,active^hit,V);
 }
 sol[h]=reference[h];
}
int main(){
 int n;cin>>n;if(n!=53)return 2;
 E marker;int old,k;cin>>old>>k>>marker.y>>marker.r[0]>>marker.r[1]>>marker.r[2];
 if(k!=1||old!=7||marker.y!=7)return 2;
 for(int i=0;i<52;i++){
  cin>>reference[i]>>k;
  if(reference[i]>=0)init|=bit(i);
  while(k--){
   E e;cin>>e.y>>e.r[0]>>e.r[1]>>e.r[2];
   bool ok=true;for(int z=0;z<3;z++)if(e.r[z]==marker.r[z])ok=false;
   if(ok)R[i].push_back(e);
  }
 }
 for(drop=0;drop<3;drop++){
  U own[3][56]{};
  for(int i=0;i<52;i++)for(const E&e:R[i])if(e.y==reference[i])
   for(int z=0;z<3;z++)if(z!=drop)own[z][e.r[z]]=bit(i);
  for(int i=0;i<52;i++)for(const E&e:R[i]){
   U x=0;for(int z=0;z<3;z++)if(z!=drop)x|=own[z][e.r[z]];
   block[i][e.y]=x;
  }
  for(lim=0;lim<=7;lim++){
   sol=reference;nodes=0;ans.clear();
   dfs(ALL^init,init,{});
   if(ans.empty())continue;
   cerr<<drop<<" "<<lim<<" "<<ans.size()<<" "<<nodes<<endl;
   for(const auto&w:ans){cout<<drop<<" "<<lim;for(int y:w)cout<<" "<<y;cout<<"\n";}
   break;
  }
  const int min[]={7,5,6},sz[]={23,1,1};
  if(lim!=min[drop]||(int)ans.size()!=sz[drop])return 3;
 }
}
