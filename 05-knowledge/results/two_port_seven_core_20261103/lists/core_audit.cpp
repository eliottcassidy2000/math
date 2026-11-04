// Independent exhaust: expand six full permutation orbits, then check every
// labelled graph against the necessary metric filters, without canonicalizing.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
using A=std::array<unsigned,7>;using U=std::uint64_t;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
A rows(unsigned code){A a{};int k=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j,++k)if((code>>k)&1){a[i]|=1u<<j;a[j]|=1u<<i;}return a;}
unsigned transform(const A&a,const std::array<int,7>&p){unsigned c=0;int k=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j,++k)if((a[p[i]]>>p[j])&1)c|=1u<<k;return c;}
bool allowed(const A&a){
 for(int i=0;i<7;++i){if(__builtin_popcount(a[i])<3)return false;for(int j=i+1;j<7;++j)if(__builtin_popcount(a[i]&a[j])>2)return false;}
 for(int center=0;center<7;++center){unsigned unseen=a[center];
  while(unseen){unsigned component=unseen&-unseen,front=component;unseen^=component;
   while(front){unsigned next=0;for(int i=0;i<7;++i)if((front>>i)&1)next|=a[i]&a[center];next&=unseen;component|=next;unseen^=next;front=next;}
   unsigned sum=0;for(int i=0;i<7;++i)if((component>>i)&1){int d=__builtin_popcount(a[i]&component);if(d>2)return false;sum+=d;}
   unsigned count=__builtin_popcount(component);
   if(sum==2*count&&count!=6)return false;
  }
 }
 return true;
}
int main(){try{
 std::map<unsigned,unsigned> owner;std::map<unsigned,unsigned> orbit_sizes;
 for(unsigned c:{120383u,120444u,120568u,122616u,378106u,380278u}){
  auto a=rows(c);need(allowed(a),"representative passes metric filters");std::array<int,7> p;std::iota(p.begin(),p.end(),0);std::set<unsigned> orbit;
  do{orbit.insert(transform(a,p));}while(std::next_permutation(p.begin(),p.end()));orbit_sizes[c]=orbit.size();
  for(auto code:orbit)need(owner.emplace(code,c).second,"representative orbits disjoint");
 }
 U mindegree=0,eligible=0,checksum=1469598103934665603ULL;
 for(unsigned code=0;code<(1u<<21);++code){auto a=rows(code);bool deg=true;for(auto v:a)deg&=__builtin_popcount(v)>=3;mindegree+=deg;
  if(allowed(a)){++eligible;need(owner.count(code)==1,"eligible graph omitted by six orbits");checksum^=code;checksum*=1099511628211ULL;}
 }
 need(eligible==owner.size(),"full orbit mass");
 std::cout<<"{\"all_graphs\":2097152,\"min_degree3_graphs\":"<<mindegree<<",\"metric_filter_graphs\":"<<eligible
 <<",\"unlabelled_representatives\":6,\"eligible_checksum_fnv64\":"<<checksum<<",\"orbit_sizes\":{";
 bool comma=false;for(auto[c,n]:orbit_sizes){if(comma)std::cout<<",";comma=true;std::cout<<"\""<<c<<"\":"<<n;}
 std::cout<<"},\"all_permutation_orbits_expanded\":true,\"eligible_set_equals_orbit_union\":true,\"status\":\"FINITE-EXACT PASS\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
