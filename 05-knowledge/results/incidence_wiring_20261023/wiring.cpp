// Full order-eight fixed-path wiring by Hamiltonian expansion, not per-tile
// canonicalization. Binary outputs are explicitly little-endian.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
void need(bool b,const string& s){if(!b)throw runtime_error(s);}
struct Node {uint32_t code;int H,aut,fibre,blue,conv;array<uint8_t,8> out;};
vector<Node> nodes;
vector<uint16_t> owner(1u<<21,65535),visits(1u<<21);
int pathv[8],current;uint64_t pathcount=0;
int cell[8][8],revcell[21];
void le(ostream& out,uint64_t x,int bytes){for(int k=0;k<bytes;k++)out.put(char(x>>(8*k)));}
void dfs(int depth,int used){
 if(depth==8){
  uint32_t mask=0;
  for(int i=0;i<8;i++)for(int j=i+2;j<8;j++)
   if(nodes[current].out[pathv[i]]>>pathv[j]&1)mask|=1u<<cell[i][j];
  need(owner[mask]==65535||owner[mask]==current,"overlapping ordinary classes");
  owner[mask]=current;visits[mask]++;pathcount++;return;
 }
 int available=(~used)&255;
 if(depth)available&=nodes[current].out[pathv[depth-1]];
 while(available){int v=__builtin_ctz(available);available&=available-1;pathv[depth]=v;dfs(depth+1,used|1<<v);}
}
array<uint8_t,8> adj(uint32_t mask){
 array<uint8_t,8>a{};for(int i=0;i<7;i++)a[i]|=1u<<(i+1);
 for(int i=0;i<8;i++)for(int j=i+2;j<8;j++)
  if(mask>>cell[i][j]&1)a[i]|=1u<<j;else a[j]|=1u<<i;
 return a;
}
uint32_t reflect(uint32_t m){uint32_t r=0;for(int i=0;i<21;i++)if(m>>i&1)r|=1u<<revcell[i];return r;}
int main(int argc,char**argv){
 need(argc==3,"usage: wiring INPUT_TEXT OUTPUT_PREFIX");
 int k=0;for(int i=0;i<8;i++)for(int j=i+2;j<8;j++)cell[i][j]=k++;
 for(int i=0;i<8;i++)for(int j=i+2;j<8;j++)revcell[cell[i][j]]=cell[7-j][7-i];
 ifstream input(argv[1]);need(bool(input),"input open");int id;Node z;
 while(input>>id>>z.code>>z.H>>z.aut>>z.fibre>>z.blue>>z.conv){
  need(id==int(nodes.size()),"class index order");z.out.fill(0);int q=27;
  for(int i=0;i<8;i++)for(int j=i+1;j<8;j++,q--)
   if(z.code>>q&1)z.out[i]|=1u<<j;else z.out[j]|=1u<<i;
  nodes.push_back(z);
 }
 need(nodes.size()==6880,"class census size");
 for(current=0;current<int(nodes.size());current++){
  auto before=pathcount;dfs(0,0);need(pathcount-before==uint64_t(nodes[current].H),"literal HP count");
 }
 vector<int>fib(6880),bluefib(6880);uint64_t blue=0,closed=0;
 map<uint64_t,uint32_t>ports,edges;
 constexpr uint32_t F=(1u<<21)-1;
 for(uint32_t m=0;m<=F;m++){
  int c=owner[m];need(c!=65535,"uncovered fixed-path presentation");
  need(visits[m]==nodes[c].aut,"path-orbit multiplicity");fib[c]++;
  uint32_t r=reflect(m);bool b=r==m;blue+=b;bluefib[c]+=b;
  need(owner[r]==nodes[c].conv,"literal converse reflection");
  need((reflect(m^F)==(r^F)),"commuting involutions");
  auto A=adj(m);int a=__builtin_popcount(A[0]),d=__builtin_popcount(A[7]);
  int target=owner[m^F];auto B=adj(m^F);
  need(__builtin_popcount(B[0])==8-a&&__builtin_popcount(B[7])==6-d,"endpoint transport");
  if(A[7]&1){
   closed++;
   uint64_t key=uint64_t(c)|(uint64_t(target)<<13)|(uint64_t(a)<<26)|(uint64_t(d)<<29)|(uint64_t(b)<<32);
   ports[key]++;
  }
  if(m<(m^F)){
   int u=min(c,nodes[c].conv),v=min(target,nodes[target].conv);if(u>v)swap(u,v);
   edges[uint64_t(u)|(uint64_t(v)<<13)|(uint64_t(b)<<26)]++;
  }
 }
 for(int c=0;c<6880;c++)need(fib[c]==nodes[c].fibre&&bluefib[c]==nodes[c].blue,"class fibre and blue inventory");
 for(auto e:edges)if(!(e.first>>26&1))need(e.second%2==0,"literal black twin multiplicity");
 string prefix=argv[2];ofstream own(prefix+".owners",ios::binary),pt(prefix+".ports",ios::binary),ed(prefix+".edges",ios::binary);
 need(bool(own)&&bool(pt)&&bool(ed),"output open");
 for(auto x:owner)le(own,x,2);
 for(auto p:ports){le(pt,p.first,8);le(pt,p.second,4);}
 for(auto e:edges){le(ed,e.first,8);le(ed,e.second,4);}
 cout<<"{\"classes\":6880,\"tilings\":"<<owner.size()<<",\"literal_paths\":"<<pathcount
 <<",\"closed_ports\":"<<closed<<",\"blue_tilings\":"<<blue
 <<",\"aggregated_port_entries\":"<<ports.size()<<",\"merged_colored_entries\":"<<edges.size()<<"}\n";
}
