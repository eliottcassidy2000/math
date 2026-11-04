#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <vector>
using U=std::uint64_t;using Mask=std::uint32_t;using Rows=std::array<unsigned,8>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
int pc(unsigned m){return __builtin_popcount(m);}
Rows decode(int n,Mask code){Rows a{};int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)if((code>>k)&1){a[i]|=1u<<j;a[j]|=1u<<i;}return a;}
Mask encode(int n,const Rows&a,const std::vector<int>&p){Mask code=0;int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)if((a[p[i]]>>p[j])&1)code|=1u<<k;return code;}
bool filters(int n,const Rows&a){
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if(pc(a[i]&a[j])>2)return false;
 for(int v=0;v<n;++v){unsigned left=a[v];while(left){unsigned todo=left&-left,component=0;bool cycle=true;while(todo){int u=__builtin_ctz(todo);todo&=todo-1;if(component&(1u<<u))continue;component|=1u<<u;unsigned near=a[u]&a[v];int degree=pc(near);if(degree>2)return false;cycle &=degree==2;todo|=near&~component;}left&=~component;if(cycle&&pc(component)!=6)return false;}}
 return true;
}
void vector_json(const std::vector<unsigned>&x){std::cout<<"[";for(size_t i=0;i<x.size();++i){if(i)std::cout<<",";std::cout<<x[i];}std::cout<<"]";}
void certificate(int n,const Rows&a){
 std::vector<int> pow3(n);int length=1;for(int i=0;i<n;++i){pow3[i]=length;length*=3;}
 std::vector<std::int64_t> value(length),next(length);value[0]=1;int m=0;
 for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)if((a[i]>>j)&1){++m;std::fill(next.begin(),next.end(),0);for(int s=0;s<length;++s)if(value[s]){if(s/pow3[i]%3<2)next[s+pow3[i]]+=value[s];if(s/pow3[j]%3<2)next[s+pow3[j]]-=value[s];}value.swap(next);}
 int first=-1;U support=0;for(int s=0;s<length;++s)if(value[s]){if(first<0)first=s;++support;}
 std::cout<<",\"polynomial_support\":"<<support<<",\"list_certificate\":";
 if(first<0){std::cout<<"null";return;}std::vector<unsigned> d;unsigned degree=0;for(int i=0;i<n;++i){d.push_back(first/pow3[i]%3);degree+=d.back();}need(degree==unsigned(m),"homogeneous certificate");std::cout<<"{\"coefficient\":"<<value[first]<<",\"exponents\":";vector_json(d);std::cout<<"}";
}
int main(int argc,char**argv){try{
 int n=argc==2?std::stoi(argv[1]):8;need(n==7||n==8,"orders7 or8");int bits=n*(n-1)/2;Mask universe=1u<<bits;std::vector<std::pair<int,int>> pairs;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)pairs.push_back({i,j});
 Rows a{};std::array<int,8> degree{};std::vector<Mask> eligible;U min3=0;Mask code=0;
 for(Mask step=0;step<universe;++step){if(step){int b=__builtin_ctz(step);auto[i,j]=pairs[b];bool adding=!((code>>b)&1);code^=1u<<b;a[i]^=1u<<j;a[j]^=1u<<i;degree[i]+=adding?1:-1;degree[j]+=adding?1:-1;}
  bool okay=true;for(int i=0;i<n;++i)okay &=degree[i]>=3;if(!okay)continue;++min3;if(filters(n,a))eligible.push_back(code);
  if(step&&step%16777216==0)std::cerr<<"step="<<step<<" survivors="<<eligible.size()<<"\n";
 }
 std::sort(eligible.begin(),eligible.end());std::unordered_set<Mask> left(eligible.begin(),eligible.end());U checksum=1469598103934665603ULL;for(Mask c:eligible){checksum^=c;checksum*=1099511628211ULL;}
 std::cout<<"{\"n\":"<<n<<",\"all_graphs\":"<<universe<<",\"min_degree3_graphs\":"<<min3<<",\"metric_filter_graphs\":"<<eligible.size()<<",\"eligible_checksum_fnv64\":"<<checksum<<",\"representatives\":[";
 U fact=1;for(int i=2;i<=n;++i)fact*=i;bool comma=false;int reps=0;
 for(Mask c:eligible){if(!left.count(c))continue;auto rows=decode(n,c);std::vector<int> p(n);std::iota(p.begin(),p.end(),0);std::unordered_set<Mask> orbit;U aut=0;
  do{Mask image=encode(n,rows,p);orbit.insert(image);aut+=image==c;}while(std::next_permutation(p.begin(),p.end()));need(orbit.size()*aut==fact,"orbit-stabilizer");for(Mask image:orbit)need(left.erase(image)==1,"disjoint complete filtered orbit");
  if(comma)std::cout<<",";comma=true;++reps;std::vector<unsigned> degrees,rv;for(int i=0;i<n;++i){degrees.push_back(pc(rows[i]));rv.push_back(rows[i]);}std::sort(degrees.begin(),degrees.end());
  std::cout<<"{\"mask\":"<<c<<",\"orbit\":"<<orbit.size()<<",\"automorphisms\":"<<aut<<",\"edges\":"<<pc(c)<<",\"rows\":";vector_json(rv);std::cout<<",\"degrees\":";vector_json(degrees);certificate(n,rows);std::cout<<"}";
 }
 need(left.empty(),"full survivor coverage");std::cout<<"],\"classes\":"<<reps<<",\"status\":\"FINITE-EXACT COMPLETE NECESSARY-FILTER UNIVERSE\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
