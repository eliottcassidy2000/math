// Complete fixed-involution codec, retaining full sidecar automorphisms.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>
using U=std::uint64_t;using W=std::uint32_t;
using Rows=std::array<unsigned,11>;using Row=std::array<W,11>;using Key=std::vector<W>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
Rows decode(int n,U c){Rows a{};int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)if((c>>k)&1)a[i]|=1u<<j;else a[j]|=1u<<i;return a;}
U encode(int n,const Rows&a,const std::vector<int>&p){U c=0;int k=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)if((a[p[i]]>>p[j])&1)c|=U(1)<<k;return c;}
std::vector<std::vector<int>> autos(int h,const Rows&a,U code){std::vector<std::vector<int>> out;std::vector<int> p(h);std::iota(p.begin(),p.end(),0);do{if(encode(h,a,p)==code)out.push_back(p);}while(std::next_permutation(p.begin(),p.end()));return out;}
unsigned transform(int h,unsigned key,const std::vector<int>&p){
 int pairs=h*(h-1)/2;unsigned rmask=key&((1u<<pairs)-1),d=key>>pairs;Rows r{};int k=0;
 for(int i=0;i<h;++i)for(int j=i+1;j<h;++j,++k)if((rmask>>k)&1){r[i]|=1u<<j;r[j]|=1u<<i;}
 unsigned out=0;k=0;for(int i=0;i<h;++i)for(int j=i+1;j<h;++j,++k)if((r[p[i]]>>p[j])&1)out|=1u<<k;
 for(int i=0;i<h;++i)if((d>>p[i])&1)out|=1u<<(pairs+i);return out;
}
Rows lift(int h,const Rows&t,unsigned key){
 Rows a{};int pairs=h*(h-1)/2;unsigned r=key&((1u<<pairs)-1),d=key>>pairs;
 auto arc=[&](int i,int j){a[i]|=1u<<j;};int k=0;
 for(int i=0;i<h;++i)for(int j=i+1;j<h;++j,++k){
  if((t[i]>>j)&1){arc(i,j);arc(h+j,h+i);}else{arc(j,i);arc(h+i,h+j);}
  if((r>>k)&1){arc(i,h+j);arc(j,h+i);}else{arc(h+j,i);arc(h+i,j);}
 }
 for(int i=0;i<h;++i){arc(i,h+i);if((d>>i)&1){arc(2*h,i);arc(h+i,2*h);}else{arc(i,2*h);arc(2*h,h+i);}}
 for(int i=0;i<2*h+1;++i)for(int j=0;j<2*h+1;++j)if(i!=j){
  int si=i==2*h?i:(i<h?i+h:i-h),sj=j==2*h?j:(j<h?j+h:j-h);
  need(((a[i]>>j)&1)==((a[sj]>>si)&1),"literal anti-involution");
 }
 return a;
}
std::vector<Row> ends(int n,const Rows&a){
 unsigned full=(1u<<n)-1;std::vector<Row> dp(full+1);Rows pred{};for(int i=0;i<n;++i)pred[i]=full^(1u<<i)^a[i];
 for(unsigned mask=1;mask<=full;++mask){unsigned vs=mask;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;unsigned rest=mask^(1u<<v),bits=rest&pred[v];W val=!rest;
  while(bits){int u=__builtin_ctz(bits);bits&=bits-1;val+=dp[rest][u];}dp[mask][v]=val;}}
 return dp;
}
struct Data{int n;Rows a;Row s{},e{};std::array<Row,11> c{};};
Data data(int n,const Rows&a){
 Data d;d.n=n;d.a=a;unsigned full=(1u<<n)-1;auto e=ends(n,a);Rows rev{};for(int i=0;i<n;++i)rev[i]=full^(1u<<i)^a[i];auto s=ends(n,rev);d.s=s[full];d.e=e[full];
 for(unsigned mask=1;mask<full;++mask){unsigned ls=mask;while(ls){int u=__builtin_ctz(ls);ls&=ls-1;if(!e[mask][u])continue;unsigned rs=full^mask;
  while(rs){int v=__builtin_ctz(rs);rs&=rs-1;d.c[u][v]+=e[mask][u]*s[full^mask][v];}}}return d;
}
std::vector<int> colors(const std::vector<Key>&keys){std::map<Key,int> m;for(auto&k:keys)m[k]=0;int next=0;for(auto&kv:m)kv.second=next++;std::vector<int> c;for(auto&k:keys)c.push_back(m.at(k));return c;}
std::vector<int> refine(const Data&d,std::vector<int> c){
 for(;;){std::vector<Key> keys;for(int u=0;u<d.n;++u){std::vector<std::tuple<int,W,W>> triples;for(int v=0;v<d.n;++v)if(v!=u)triples.emplace_back(c[v],d.c[u][v],d.c[v][u]);
  std::sort(triples.begin(),triples.end());Key k{W(c[u])};for(auto[col,a,b]:triples){k.push_back(col);k.push_back(a);k.push_back(b);}keys.push_back(k);}
  auto nxt=colors(keys);if(nxt==c)return c;c=nxt;
 }
}
struct Canon{Key key;U lo=0,hi=0,aut=0,leaves=0,nodes=0;};
Canon canonical(const Data&d){
 std::vector<Key> start;for(int i=0;i<d.n;++i)start.push_back({d.s[i],d.e[i]});Canon out;
 std::function<void(std::vector<int>)> visit=[&](std::vector<int> c){
  ++out.nodes;c=refine(d,c);int cells=*std::max_element(c.begin(),c.end())+1;std::vector<std::vector<int>> groups(cells);for(int i=0;i<d.n;++i)groups[c[i]].push_back(i);
  int cell=-1;for(int i=0;i<cells;++i)if(groups[i].size()>1){cell=i;break;}
  if(cell<0){std::vector<int> p;for(auto&g:groups)p.push_back(g[0]);Key key;for(int u:p){key.push_back(d.s[u]);key.push_back(d.e[u]);}for(int u:p)for(int v:p)key.push_back(d.c[u][v]);++out.leaves;
   if(out.aut&&key>out.key)return;U code=encode(d.n,d.a,p);
   if(!out.aut||key<out.key){out.key=key;out.lo=out.hi=code;out.aut=1;}else{out.lo=std::min(out.lo,code);out.hi=std::max(out.hi,code);++out.aut;}return;
  }
  for(int v:groups[cell]){auto split=c;for(int u=0;u<d.n;++u)split[u]=2*c[u]+(c[u]==cell&&u!=v);visit(split);}
 };
 visit(colors(start));return out;
}
void histogram(const std::map<U,U>&h){std::cout<<"{";bool comma=false;for(auto[k,v]:h){if(comma)std::cout<<",";comma=true;std::cout<<"\""<<k<<"\":"<<v;}std::cout<<"}";}
int main(int argc,char**argv){try{
 need(argc==2,"usage: census core-banks.txt");std::ifstream input(argv[1]);need(bool(input),"core banks");int maxh;input>>maxh;need(maxh>=1&&maxh<=5,"half-order bound");
 std::cout<<"{\"orders\":[";
 for(int h=1;h<=maxh;++h){int count;input>>count;std::map<Key,U> seen;std::map<U,U> center_hist,full_hist,index_hist;U reps=0,source_mass=0,classbad=0,autbad=0,leaves=0,nodes=0,maxleaf=0,balanced=0;
  std::vector<std::array<U,3>> collisions;U fact=1;for(int i=2;i<=h;++i)fact*=i;
  for(int ti=0;ti<count;++ti){U corecode,expected_aut;input>>corecode>>expected_aut;need(bool(input),"complete core record");auto t=decode(h,corecode);auto group=autos(h,t,corecode);need(group.size()==expected_aut,"literal core automorphisms");
   for(unsigned key=0;key<(1u<<(h*(h-1)/2+h));++key){bool keep=true;U stabilizer=0;
    for(auto&p:group){auto image=transform(h,key,p);if(image<key){keep=false;break;}stabilizer+=image==key;}if(!keep)continue;
    ++reps;++center_hist[stabilizer];need(fact%stabilizer==0,"source orbit size");source_mass+=fact/stabilizer;
    auto d=data(2*h+1,lift(h,t,key));auto result=canonical(d);balanced+=d.s==d.e;leaves+=result.leaves;nodes+=result.nodes;maxleaf=std::max(maxleaf,result.leaves);
    ++full_hist[result.aut];need(result.aut%stabilizer==0,"centralizer divides full group");++index_hist[result.aut/stabilizer];
    if(result.lo!=result.hi){++autbad;collisions.push_back({1,result.lo,result.hi});}
    auto ins=seen.emplace(result.key,result.lo);if(!ins.second){++classbad;collisions.push_back({2,ins.first->second,result.lo});}
    if(reps%25000==0)std::cerr<<"half="<<h<<" representatives="<<reps<<"\n";
   }
  }
  need(source_mass==(U(1)<<(h*h)),"source labelled mass");
  if(h>1)std::cout<<",";std::cout<<"{\"n\":"<<2*h+1<<",\"classes\":"<<reps<<",\"source_labelled_mass\":"<<source_mass<<",\"sidecar_classes\":"<<seen.size()
   <<",\"class_collisions\":"<<classbad<<",\"extra_sidecar_symmetries\":"<<autbad<<",\"endpoint_balanced_classes\":"<<balanced
   <<",\"canonical_leaves\":"<<leaves<<",\"canonical_nodes\":"<<nodes<<",\"max_leaves\":"<<maxleaf<<",\"centralizer_histogram\":";histogram(center_hist);
  std::cout<<",\"full_automorphism_histogram\":";histogram(full_hist);std::cout<<",\"anti_involution_histogram\":";histogram(index_hist);
  std::cout<<",\"collisions\":[";bool comma=false;for(auto x:collisions){if(comma)std::cout<<",";comma=true;std::cout<<"["<<x[0]<<","<<x[1]<<","<<x[2]<<"]";}std::cout<<"]}";
 }
 std::cout<<"],\"scope\":\"All self-converse tournaments at odd orders3..11 (as requested max half-order); not all tournaments of order11\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
