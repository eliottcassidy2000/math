// Search the unchanged physical boundary of the reference patch.
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
int main(int argc,char **argv){
 int maxSize=argc==2?stoi(argv[1]):12;
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
 for(drop=3;drop<4;drop++){
  U own[3][56]{};
  for(int i=0;i<52;i++)for(const E&e:R[i])if(e.y==reference[i])
   for(int z=0;z<3;z++)if(z!=drop)own[z][e.r[z]]=bit(i);
  for(int i=0;i<52;i++)for(const E&e:R[i]){
   U x=0;for(int z=0;z<3;z++)if(z!=drop)x|=own[z][e.r[z]];
   block[i][e.y]=x;
  }

  for(int h=0;h<52;h++) {
    vector<E> next;
    for(const E& e:R[h]){
      bool ok=(reference[h]>=0&&reference[h]!=e.y);
      for(int z=0;z<3;z++)if(!own[z][e.r[z]])ok=false;
      if(ok)next.push_back(e);
    }
    R[h]=next;
  }
  for(lim=4;lim<=maxSize;lim++){
   ans.clear();nodes=0;
   sol=reference;
   for(int h=0;h<52;h++)for(const E&e:R[h]){
    U hit=block[h][e.y]|bit(h);
    if(__builtin_popcountll(hit)>lim)continue;
    sol[h]=e.y;
    array<U,3> V{};
    for(int z=0;z<3;z++)V[z]|=bit(e.r[z]);
    dfs(hit^bit(h),init^hit,V);
    sol[h]=reference[h];
   }
   cerr<<lim<<" "<<ans.size()<<" "<<nodes<<endl;
   if(ans.empty())continue;
   for(const auto&w:ans){cout<<lim;for(int y:w)cout<<" "<<y;cout<<"\n";}
   return 0;
  }
 }
 cout<<"No same-boundary trade through size "<<maxSize<<endl;
}
