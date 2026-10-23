#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
// Discovery only: q=3, m=5, complete 2^18 marked Hamilton-cycle frames.
constexpr int N=15,Q=3,M=5;
struct Orbit{std::vector<std::pair<int,int>> edges;int forced=-1;};
std::vector<Orbit> orbits;
std::array<uint32_t,N> out;
std::array<int,N> path;
using Key=std::array<uint64_t,2>;
std::set<Key> keys;
std::vector<uint64_t> witnesses;
uint64_t nodes=0, words=0;
void demand(bool x,const char*s){if(!x)throw std::runtime_error(s);}
bool arc(int x,int y){return (out[x]>>y)&1;}
void initialize(){
 bool seen[N][N]{};
 for(int a=0;a<N;++a)for(int b=a+1;b<N;++b)if(!seen[a][b]){
  Orbit o;
  for(int k=0;k<M;++k){int x=(a+Q*k)%N,y=(b+Q*k)%N;
   seen[std::min(x,y)][std::max(x,y)]=true;o.edges.push_back({x,y});
   int f=(x+1)%N==y?1:((y+1)%N==x?0:-1);
   if(f>=0){demand(o.forced<0||o.forced==f,"cycle contradiction");o.forced=f;}}
  orbits.push_back(o);
 }
 demand(orbits.size()==21,"21 orbits");
 demand(std::count_if(orbits.begin(),orbits.end(),[](auto&o){return o.forced<0;})==18,"18 free bits");
}
void frame(unsigned mask){out.fill(0);int k=0;
 for(auto&o:orbits){int bit=o.forced<0?int(mask>>k++&1):o.forced;
  for(auto[x,y]:o.edges){if(!bit)std::swap(x,y);out[x]|=1u<<y;}}
}
bool symmetric(){std::array<int,N> hist{};for(auto r:out)++hist[std::popcount(r)];
 for(int i=0;i<N;++i)if(hist[i]!=hist[N-1-i])return false;return true;}
void record(){++words;Key key{};bool blue=true;int pos=0;uint64_t encoded=0;
 for(int i=0;i<N;++i){encoded|=uint64_t(path[i])<<(4*i);
  for(int j=i+1;j<N;++j,++pos){bool a=arc(path[i],path[j]);
   if(a)key[pos/64]|=uint64_t(1)<<(pos%64);
   if(a==arc(path[N-1-i],path[N-1-j]))blue=false;}}
 if(!blue&&keys.insert(key).second)witnesses.push_back(encoded);
}
void dfs(int at,unsigned used,int depth){
 if(witnesses.size()==7)return;
 if(++nodes>10000000)throw std::runtime_error("discovery cap exhausted");
 if(depth==N){record();return;}
 unsigned next=out[at]&~used;
 while(next&&witnesses.size()<7){unsigned bit=next&-next;next-=bit;int j=std::countr_zero(bit);
  path[depth]=j;dfs(j,used|bit,depth+1);}
}
int main(int argc,char**argv){
 demand(argc==2,"raw output path");initialize();std::ofstream dst(argv[1],std::ios::binary);
 demand(bool(dst),"output open");dst.write("D15Q302!",8);
 uint64_t total=0,maxnodes=0,maxwords=0;unsigned rejected=0,accepted=0;
 for(unsigned mask=0;mask<(1u<<18);++mask){frame(mask);
  if(!symmetric()){dst.put(0);++rejected;continue;}
  keys.clear();witnesses.clear();nodes=words=0;path[0]=0;dfs(0,1,1);
  demand(witnesses.size()==7,"missing seven words");dst.put(1);
  __uint128_t source=0;int edge=0;
  for(int a=0;a<N;++a)for(int b=a+1;b<N;++b,++edge)
   if(arc(a,b))source|=__uint128_t(1)<<edge;
  for(uint64_t w:witnesses){
   for(int k=0;k<8;++k)dst.put(char(w>>(8*k)&255));
   __uint128_t kept=0;
   for(int k=0;k<N-1;++k){int a=w>>(4*k)&15,b=w>>(4*(k+1))&15;
    if(a>b)std::swap(a,b);int e=a*(2*N-a-1)/2+b-a-1;
    kept|=__uint128_t(1)<<e;}
   __uint128_t target=source^((__uint128_t(1)<<105)-1)^kept;
   for(int k=0;k<14;++k)dst.put(char(target>>(8*k)&255));
  }
  ++accepted;total+=nodes;maxnodes=std::max(maxnodes,nodes);maxwords=std::max(maxwords,words);
 }
 demand(bool(dst),"output complete");
 std::cout<<"frames262144 rejected "<<rejected<<" certified "<<accepted<<" nodes "<<total<<" maxnodes "<<maxnodes<<" maxwords "<<maxwords<<'\n';
}
