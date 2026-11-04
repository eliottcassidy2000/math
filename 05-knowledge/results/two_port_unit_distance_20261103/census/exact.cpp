#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>
using U=std::uint64_t;
using Key=std::vector<U>;
using Row=std::array<U,9>;

struct Data {
 int n; std::array<unsigned,9> a{}; Row s{},e{}; std::array<Row,9> c{};
};
std::vector<Row> ends(int n,const std::array<unsigned,9>& a) {
 std::vector<Row> dp(1<<n);
 std::array<unsigned,9> pred{};
 for(int u=0;u<n;++u)for(int v=0;v<n;++v)if((a[u]>>v)&1)pred[v]|=1u<<u;
 for(unsigned mask=1;mask<(1u<<n);++mask) for(int v=0;v<n;++v)if(mask&(1u<<v)) {
  unsigned rest=mask^(1u<<v), bits=rest&pred[v];
  if(!rest)dp[mask][v]=1;
  while(bits){int u=__builtin_ctz(bits);bits&=bits-1;dp[mask][v]+=dp[rest][u];}
 }
 return dp;
}
Data make(int n,U code){
 Data d;d.n=n;int bit=0;
 for(int u=0;u<n;++u)for(int v=u+1;v<n;++v,++bit)
  if((code>>bit)&1)d.a[u]|=1u<<v;else d.a[v]|=1u<<u;
 auto e=ends(n,d.a);std::array<unsigned,9> rev{};unsigned full=(1u<<n)-1;
 for(int u=0;u<n;++u)rev[u]=full^(1u<<u)^d.a[u];
 auto s=ends(n,rev);d.e=e[full];d.s=s[full];
 for(unsigned mask=1;mask<full;++mask)
  for(int u=0;u<n;++u)if(e[mask][u])
   for(int v=0;v<n;++v)if(s[full^mask][v])d.c[u][v]+=e[mask][u]*s[full^mask][v];
 return d;
}
std::vector<int> colors_of(const std::vector<Key>& keys){
 std::map<Key,int> ids;for(auto &k:keys)ids[k]=0;
 int id=0;for(auto &p:ids)p.second=id++;
 std::vector<int> out;for(auto &k:keys)out.push_back(ids.at(k));return out;
}
struct Canon {Key key; U mincode=0,maxcode=0,aut=0,leaves=0;};
Canon canon(const Data& d){
 int n=d.n;std::vector<Key> keys;
 for(int u=0;u<n;++u)keys.push_back({d.s[u],d.e[u]});
 auto colors=colors_of(keys);
 for(;;){
  keys.clear();
  for(int u=0;u<n;++u){
   std::vector<std::tuple<int,U,U>> triples;
   for(int v=0;v<n;++v)if(u!=v)triples.emplace_back(colors[v],d.c[u][v],d.c[v][u]);
   std::sort(triples.begin(),triples.end());Key k{U(colors[u])};
   for(auto [col,a,b]:triples){k.push_back(col);k.push_back(a);k.push_back(b);}keys.push_back(k);
  }
  auto next=colors_of(keys);if(next==colors)break;colors=next;
 }
 std::vector<std::vector<int>> cells(*std::max_element(colors.begin(),colors.end())+1);
 for(int u=0;u<n;++u)cells[colors[u]].push_back(u);
 std::vector<int> order;Canon result;
 std::function<void(int)> visit=[&](int at){
  if(at==(int)cells.size()){
   Key key;for(int u:order){key.push_back(d.s[u]);key.push_back(d.e[u]);}
   for(int u:order)for(int v:order)key.push_back(d.c[u][v]);++result.leaves;
   if(result.aut && key>result.key)return;
   U code=0;int bit=0;
   for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)if((d.a[order[i]]>>order[j])&1)code|=U(1)<<bit;
   if(!result.aut || key<result.key){result.key=key;result.mincode=result.maxcode=code;result.aut=1;}
   else{result.mincode=std::min(result.mincode,code);result.maxcode=std::max(result.maxcode,code);++result.aut;}
   return;
  }
  auto cell=cells[at];
  do {auto old=order.size();order.insert(order.end(),cell.begin(),cell.end());visit(at+1);order.resize(old);}
  while(std::next_permutation(cell.begin(),cell.end()));
 };
 visit(0);return result;
}
int main(int argc,char**argv){
 try{
  if(argc!=2)throw std::runtime_error("expected temporary input path");
  std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open input");
  int orders;input>>orders;U collisions=0;
  std::cout<<"{\"orders\":[";
  for(int at=0;at<orders;++at){
   int n;U count;input>>n>>count;std::map<Key,U> seen;U aut_bad=0,total=0,maximum=0;
   for(U index=0;index<count;++index){
    U code,expected_aut;input>>code>>expected_aut;
    if(!input)throw std::runtime_error("truncated input");
    auto result=canon(make(n,code));total+=result.leaves;maximum=std::max(maximum,result.leaves);
    if(result.aut!=expected_aut)++aut_bad;
    if(result.mincode!=result.maxcode){++collisions;std::cerr<<"AUT COLLISION "<<n<<" "<<result.mincode<<" "<<result.maxcode<<"\n";}
    auto inserted=seen.emplace(result.key,result.mincode);
    if(!inserted.second){++collisions;std::cerr<<"CLASS COLLISION "<<n<<" "<<inserted.first->second<<" "<<result.mincode<<"\n";}
   }
   if(at)std::cout<<",";
   std::cout<<"{\"n\":"<<n<<",\"classes\":"<<count<<",\"sidecar_classes\":"<<seen.size()
    <<",\"automorphism_group_disagreements\":"<<aut_bad<<",\"canonical_leaves\":"<<total<<",\"max_leaves\":"<<maximum<<"}";
   std::cerr<<"completed order "<<n<<": "<<count<<" classes\n";
  }
  std::cout<<"],\"collision_count\":"<<collisions<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
