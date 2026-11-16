#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>
using A=std::array<int,9>;
unsigned long long key(const A&a){
 A p;std::iota(p.begin(),p.end(),0);
 auto deg=[&](int i){return __builtin_popcount(a[i]);};
 std::sort(p.begin(),p.end(),[&](int i,int j){return deg(i)==deg(j)?i<j:deg(i)<deg(j);});
 unsigned long long best=~0ULL;
 auto rec=[&](auto&&self,int i)->void{
  if(i==9){
   unsigned long long k=0;
   for(int t=0;t<9;t++)for(int q=t+1;q<9;q++)k=(k<<1)|((a[p[t]]>>p[q])&1);
   best=std::min(best,k);return;
  }
  int j=i+1;while(j<9&&deg(p[j])==deg(p[i]))j++;
  do{self(self,j);}while(std::next_permutation(p.begin()+i,p.begin()+j));
 };
 rec(rec,0);return best;
}
int main(){
 std::unordered_map<unsigned long long,std::array<int,3>> targets;
 const std::array<std::array<int,8>,2> hosts={{{2,16,3,23,165,79,31,111},{16,5,17,7,74,31,175,63}}};
 for(int hostidx=0;hostidx<2;hostidx++)for(int mask=0;mask<256;mask++){
  auto h=hosts[hostidx];
  std::vector<int>s;for(int i=0;i<8;i++)if(mask>>i&1)s.push_back(i);
  int k=s.size();for(int p=0;p+k<=9;p++){
   A a{}; for(int i=0;i<9;i++)for(int j=i+1;j<9;j++){
    if(i>=p&&j<p+k&&!(h[s[i-p]]>>s[j-p]&1))a[j]|=1<<i;else a[i]|=1<<j;
   }
   auto z=key(a); auto it=targets.find(z);
   if(it==targets.end() || __builtin_popcount(it->second[1])<k)targets[z]={hostidx,mask,p};
  }
 }
 std::cout<<"allowed patterns "<<targets.size()<<"\n";
 A a{};std::vector<std::pair<int,int>> edges;
 for(int i=0;i<9;i++)for(int j=i+1;j<9;j++){edges.emplace_back(i,j);if(i/3==j/3&&i%3==0&&j%3==2)a[j]|=1<<i;else a[i]|=1<<j;}
 for(int k=0;k<=5;k++){
  int hit=0,bestk=0; std::vector<std::pair<int,int>> flips;
  auto rec=[&](auto&&self,int l,int e)->void{
   if(l==0){
    auto it=targets.find(key(a));
    if(it!=targets.end()){
     auto [ii,mm,pp]=it->second;
     if(!hit || __builtin_popcount(mm)>bestk){
      bestk=__builtin_popcount(mm);
      for(auto [u,v]:flips)std::cout<<u<<"-"<<v<<" ";
      auto [i,mask,p]=it->second;std::cout<<" host "<<i<<" mask "<<mask<<" before "<<p<<"\n";}
     hit++;
    }
    return;
   }
   for(int t=e;t+l<=int(edges.size());t++){
    auto [i,j]=edges[t];a[i]^=1<<j;a[j]^=1<<i;flips.emplace_back(i,j);self(self,l-1,t+1);flips.pop_back();a[i]^=1<<j;a[j]^=1<<i;
   }
  };
  rec(rec,k,0);std::cout<<k<<" flips: "<<hit<<" hits\n"<<std::flush;if(hit)break;
 }
}
