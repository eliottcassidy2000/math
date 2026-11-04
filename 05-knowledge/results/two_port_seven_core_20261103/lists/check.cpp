// Independent complete 3-list census from a fixed five-colour palette.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
using A=std::array<unsigned,7>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
A graph(unsigned code){A a{};int k=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j,++k)if((code>>k)&1){a[i]|=1u<<j;a[j]|=1u<<i;}return a;}
bool colour(const A&a,const A&lists,unsigned done,A&chosen){
 if(done==127)return true;
 int pick=-1,best=6;unsigned allowed=0;
 for(int u=0;u<7;++u)if(!((done>>u)&1)){
  unsigned colors=lists[u],ns=a[u]&done;
  while(ns){int v=__builtin_ctz(ns);ns&=ns-1;colors&=~chosen[v];}
  int size=__builtin_popcount(colors);if(!size)return false;
  if(size<best){best=size;pick=u;allowed=colors;}
 }
 while(allowed){unsigned c=allowed&-allowed;allowed-=c;chosen[pick]=c;if(colour(a,lists,done|(1u<<pick),chosen))return true;}
 chosen[pick]=0;return false;
}
int main(){try{
 std::vector<unsigned> choices;for(unsigned m=0;m<32;++m)if(__builtin_popcount(m)==3)choices.push_back(m);
 std::cout<<"{\"palette\":5,\"list_size\":3,\"normalized_vertex\":0,\"normalized_list\":7,\"graphs\":[";
 bool comma=false;
 for(unsigned code:{120383u,120444u,120568u,378106u}){
  auto a=graph(code);U bad=0;std::vector<A> failures;U witness=1469598103934665603ULL;
  for(unsigned index=0;index<1000000;++index){unsigned x=index;A lists{},chosen{};lists[0]=7;
   for(int u=1;u<7;++u){lists[u]=choices[x%10];x/=10;}
   bool ok=colour(a,lists,0,chosen);
   if(!ok){++bad;failures.push_back(lists);}
   else{for(int u=0;u<7;++u){need((lists[u]&chosen[u])!=0&&__builtin_popcount(chosen[u])==1,"list witness");
      for(int v=0;v<u;++v)need(!((a[u]>>v)&1)||chosen[u]!=chosen[v],"edge witness");
      witness^=chosen[u];witness*=1099511628211ULL;}}
  }
  if(comma)std::cout<<",";comma=true;std::cout<<"{\"code\":"<<code<<",\"degrees\":[";
  for(int u=0;u<7;++u){if(u)std::cout<<",";std::cout<<__builtin_popcount(a[u]);}
  std::cout<<"],\"normalized_lists\":1000000,\"represented_lists\":10000000,\"bad\":"<<bad<<",\"witness_checksum_fnv64\":"<<witness<<",\"bad_lists\":[";
  for(std::size_t j=0;j<failures.size();++j){if(j)std::cout<<",";std::cout<<"[";for(int u=0;u<7;++u){if(u)std::cout<<",";std::cout<<failures[j][u];}std::cout<<"]";}
  std::cout<<"]}";
 }
 std::cout<<"],\"status\":\"FINITE-EXACT\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
