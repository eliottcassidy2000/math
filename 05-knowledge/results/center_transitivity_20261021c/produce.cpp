#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
// Discovery only. check.py reconstructs every source and consumes literal words.
using Key=std::array<uint64_t,5>;
std::array<uint32_t,25> out;
std::array<int,25> path;
std::set<Key> keys;
std::vector<std::array<int,25>> witnesses;
uint64_t nodes, words;
int diff(int x,int y,int kind){
  if(kind==0)return (y-x+25)%25;
  return ((y%5-x%5+5)%5)+5*((y/5-x/5+5)%5);
}
bool arc(int x,int y){return (out[x]>>y)&1;}
void record(){
  ++words;
  Key key{}; bool blue=true;int pos=0;
  for(int i=0;i<25;++i)for(int j=i+1;j<25;++j,++pos){
    bool a=arc(path[i],path[j]);
    if(a)key[pos/64]|=uint64_t(1)<<(pos%64);
    if(a==arc(path[24-i],path[24-j]))blue=false;
  }
  if(!blue && keys.insert(key).second)witnesses.push_back(path);
}
void dfs(int at,uint32_t used,int depth){
  if(witnesses.size()==7)return;
  if(++nodes>10000000)throw std::runtime_error("discovery cap exhausted");
  if(depth==25){record();return;}
  uint32_t next=out[at]&~used;
  while(next && witnesses.size()<7){
    int j=std::countr_zero(next);next&=next-1;
    path[depth]=j;dfs(j,used|(uint32_t(1)<<j),depth+1);
  }
}
int main(int argc,char**argv){
  if(argc!=2)throw std::runtime_error("output path required");
  std::ofstream dst(argv[1]);
  uint64_t total=0,maxnodes=0,maxwords=0;
  for(int kind=0;kind<2;++kind){
    std::vector<int> reps;
    for(int d=1;d<25;++d)if(d<diff(d,0,kind))reps.push_back(d);
    if(reps.size()!=12)throw std::runtime_error("bad inverse pairs");
    for(int mask=0;mask<4096;++mask){
      uint32_t connection=0;
      for(int i=0;i<12;++i){int d=reps[i];if((mask>>i)&1)d=diff(d,0,kind);connection|=uint32_t(1)<<d;}
      out.fill(0);
      for(int x=0;x<25;++x)for(int y=0;y<25;++y)if((connection>>diff(x,y,kind))&1)out[x]|=uint32_t(1)<<y;
      keys.clear();witnesses.clear();nodes=words=0;path[0]=0;dfs(0,1,1);
      if(witnesses.size()!=7)throw std::runtime_error("missing seven witnesses");
      dst<<kind<<' '<<mask;
      for(const auto&p:witnesses){dst<<' ';for(int v:p)dst<<char('A'+v);}
      dst<<' '<<nodes<<' '<<words<<'\n';
      total+=nodes;maxnodes=std::max(maxnodes,nodes);maxwords=std::max(maxwords,words);
    }
  }
  std::cout<<"frames8192 nodes "<<total<<" maxnodes "<<maxnodes<<" maxwords "<<maxwords<<'\n';
}
