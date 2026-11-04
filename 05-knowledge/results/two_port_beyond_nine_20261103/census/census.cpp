// Exact two-stage census. Endpoint keys are stored in full, never as hashes.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using U=std::uint64_t;
using W=std::uint32_t;
constexpr int N=10, FULL=(1<<N)-1;
using Row=std::array<W,N>;
using Rows=std::array<unsigned,N>;
using Table=std::array<Row,1<<N>;
using Key=std::vector<W>;
struct Record { std::array<W,2*N> key{}; W address=0; };
struct Data {Rows a{};Row s{},e{};std::array<Row,N> c{};};
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}

std::vector<U> read64(const std::string&path){
 std::ifstream f(path,std::ios::binary);need(bool(f),"open parents");std::vector<U> out;
 for(;;){int c=f.get();if(c==EOF)break;U x=unsigned(c);
  for(int j=1;j<8;++j){c=f.get();need(c!=EOF,"partial parent");x|=U(unsigned(c))<<(8*j);}out.push_back(x);}
 return out;
}
Rows decode(U code,unsigned cut){
 Rows a{};int bit=0;
 for(int i=0;i<9;++i)for(int j=i+1;j<9;++j,++bit)
  if((code>>bit)&1)a[i]|=1u<<j;else a[j]|=1u<<i;
 for(int i=0;i<9;++i)if((cut>>i)&1)a[9]|=1u<<i;else a[i]|=1u<<9;
 return a;
}
Table ends(const Rows&a){
 Table t{};Rows pred{};
 for(int v=0;v<N;++v)pred[v]=FULL^(1u<<v)^a[v];
 for(unsigned mask=1;mask<=FULL;++mask){
  unsigned vertices=mask;
  while(vertices){int v=__builtin_ctz(vertices);vertices&=vertices-1;
   unsigned rest=mask^(1u<<v), bits=rest&pred[v];W value=!rest;
   while(bits){int u=__builtin_ctz(bits);bits&=bits-1;value+=t[rest][u];}
   t[mask][v]=value;
  }
 }
 return t;
}
Data make(const Rows&a,bool covers){
 Data d;d.a=a;auto e=ends(a);Rows rev{};
 for(int u=0;u<N;++u)rev[u]=FULL^(1u<<u)^a[u];
 auto s=ends(rev);d.e=e[FULL];d.s=s[FULL];
 if(covers)for(unsigned mask=1;mask<FULL;++mask){
  unsigned left=mask;
  while(left){int u=__builtin_ctz(left);left&=left-1;W eu=e[mask][u];if(!eu)continue;
   unsigned right=FULL^mask;
   while(right){int v=__builtin_ctz(right);right&=right-1;d.c[u][v]+=eu*s[FULL^mask][v];}
  }
 }
 return d;
}
U code_at(const Rows&a,const std::vector<int>&p){
 U code=0;int bit=0;
 for(int i=0;i<N;++i)for(int j=i+1;j<N;++j,++bit)if((a[p[i]]>>p[j])&1)code|=U(1)<<bit;
 return code;
}
std::vector<int> colors_of(const std::vector<Key>&keys){
 std::map<Key,int> ids;for(auto&k:keys)ids[k]=0;int id=0;for(auto&p:ids)p.second=id++;
 std::vector<int> out;for(auto&k:keys)out.push_back(ids.at(k));return out;
}
struct Canon {Key key;U mincode=0,maxcode=0,aut=0,leaves=0;};
Canon canon(const Data&d){
 std::vector<Key> keys;for(int u=0;u<N;++u)keys.push_back({d.s[u],d.e[u]});
 auto colors=colors_of(keys);
 for(;;){
  keys.clear();for(int u=0;u<N;++u){
   std::vector<std::tuple<int,W,W>> triples;
   for(int v=0;v<N;++v)if(v!=u)triples.emplace_back(colors[v],d.c[u][v],d.c[v][u]);
   std::sort(triples.begin(),triples.end());Key k{W(colors[u])};
   for(auto[col,a,b]:triples){k.push_back(col);k.push_back(a);k.push_back(b);}keys.push_back(k);
  }
  auto next=colors_of(keys);if(next==colors)break;colors=next;
 }
 std::vector<std::vector<int>> cells(*std::max_element(colors.begin(),colors.end())+1);
 for(int u=0;u<N;++u)cells[colors[u]].push_back(u);
 std::vector<int> order;Canon out;
 std::function<void(int)> walk=[&](int at){
  if(at==int(cells.size())){
   Key key;for(int u:order){key.push_back(d.s[u]);key.push_back(d.e[u]);}
   for(int u:order)for(int v:order)key.push_back(d.c[u][v]);++out.leaves;
   if(out.aut&&key>out.key)return;
   U code=code_at(d.a,order);
   if(!out.aut||key<out.key){out.key=key;out.mincode=out.maxcode=code;out.aut=1;}
   else{out.mincode=std::min(out.mincode,code);out.maxcode=std::max(out.maxcode,code);++out.aut;}
   return;
  }
  auto cell=cells[at];do{auto old=order.size();order.insert(order.end(),cell.begin(),cell.end());walk(at+1);order.resize(old);}
  while(std::next_permutation(cell.begin(),cell.end()));
 };
 walk(0);return out;
}

int main(int argc,char**argv){try{
 need(argc==3||argc==4,"usage: census parents.u64 flags.varint [limit]");
 U limit=argc==4?std::stoull(argv[3]):9733056;
 auto parents=read64(argv[1]);need(parents.size()==191536,"parent count");
 std::ifstream f(argv[2],std::ios::binary);need(bool(f),"open flags");
 std::vector<Record> records;records.reserve(limit);
 U address=0,previous=0,balanced=0,endpoint_ties=0;
 for(;;){int c=f.get();if(c==EOF)break;U delta=0;int shift=0;
  for(;;){need(shift<63,"varint overflow");delta|=U(c&127)<<shift;if(!(c&128))break;shift+=7;c=f.get();need(c!=EOF,"partial varint");}
  address+=delta;need(records.empty()||address>previous,"addresses not increasing");previous=address;
  need((address>>9)<parents.size()&&address<(U(1)<<27),"address bounds");
  auto d=make(decode(parents[address>>9],address&511),false);
  balanced+=d.s==d.e;
  std::array<std::pair<W,W>,N> pairs;for(int u=0;u<N;++u)pairs[u]={d.s[u],d.e[u]};
  std::sort(pairs.begin(),pairs.end());endpoint_ties+=std::adjacent_find(pairs.begin(),pairs.end())!=pairs.end();
  Record r;r.address=address;for(int u=0;u<N;++u){r.key[2*u]=pairs[u].first;r.key[2*u+1]=pairs[u].second;}records.push_back(r);
  if(records.size()%1000000==0)std::cerr<<"endpoint_classes="<<records.size()<<"\n";
  if(records.size()==limit)break;
 }
 need(records.size()==limit,"short bank");
 if(argc==3)need(f.peek()==EOF,"extra bank entries");
 std::sort(records.begin(),records.end(),[](const Record&a,const Record&b){return a.key!=b.key?a.key<b.key:a.address<b.address;});
 U groups=0,nontrivial_groups=0,maxgroup=0,covered=0,leaves=0,maxleaves=0,class_collision=0,aut_collision=0,mass=0;
 std::map<U,U> aut_hist;std::vector<std::array<U,3>> collisions;
 for(std::size_t i=0;i<records.size();){std::size_t j=i+1;while(j<records.size()&&records[j].key==records[i].key)++j;
  ++groups;maxgroup=std::max(maxgroup,U(j-i));nontrivial_groups+=j-i>1;
  bool ties=false;for(int u=1;u<N;++u)ties|=records[i].key[2*u]==records[i].key[2*u-2]&&records[i].key[2*u+1]==records[i].key[2*u-1];
  if(j==i+1&&!ties){++aut_hist[1];mass+=3628800;i=j;continue;}
  std::map<Key,U> seen;
  for(std::size_t k=i;k<j;++k){U a=records[k].address;auto d=make(decode(parents[a>>9],a&511),true);auto out=canon(d);
   ++covered;leaves+=out.leaves;maxleaves=std::max(maxleaves,out.leaves);++aut_hist[out.aut];need(3628800%out.aut==0,"group order");mass+=3628800/out.aut;
   if(out.mincode!=out.maxcode){++aut_collision;collisions.push_back({1,out.mincode,out.maxcode});}
   auto ins=seen.emplace(out.key,out.mincode);if(!ins.second){++class_collision;collisions.push_back({2,ins.first->second,out.mincode});}
  }
  i=j;
 }
 std::cout<<"{\"order\":10,\"classes\":"<<records.size()<<",\"complete_bank\":"<<(argc==3?"true":"false")
 <<",\"endpoint_balanced_classes\":"<<balanced<<",\"classes_with_endpoint_ties\":"<<endpoint_ties
 <<",\"endpoint_multiset_groups\":"<<groups<<",\"nontrivial_endpoint_groups\":"<<nontrivial_groups<<",\"max_endpoint_group\":"<<maxgroup
 <<",\"full_cover_canonicalizations\":"<<covered<<",\"cover_canonical_leaves\":"<<leaves<<",\"max_cover_leaves\":"<<maxleaves
 <<",\"sidecar_class_collisions\":"<<class_collision<<",\"sidecar_automorphism_collisions\":"<<aut_collision
 <<",\"sidecar_labelled_mass\":"<<mass<<",\"sidecar_automorphism_histogram\":{";
 bool comma=false;for(auto[g,num]:aut_hist){if(comma)std::cout<<",";comma=true;std::cout<<"\""<<g<<"\":"<<num;}
 std::cout<<"},\"collisions\":[";comma=false;for(auto r:collisions){if(comma)std::cout<<",";comma=true;std::cout<<"["<<r[0]<<","<<r[1]<<","<<r[2]<<"]";}
 std::cout<<"]}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
