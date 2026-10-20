#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
struct Action{uint32_t flip,sign,tab[4][256]{};};
int main(int argc,char**argv){
 if(argc!=4)throw runtime_error("input output mode");ifstream in(argv[1]);ofstream out(argv[2]);int d,k;in>>d>>k;vector<Action>a(k);
 for(auto&r:a){in>>r.flip>>r.sign;vector<unsigned>im(d);for(auto&v:im)in>>v;for(int b=0;b<4;b++)for(unsigned x=1;x<256;x++){unsigned i=8*b+__builtin_ctz(x);r.tab[b][x]=r.tab[b][x&(x-1)]|(i<(unsigned)d?1u<<im[i]:0);}}
 bool T=string(argv[3])=="T";vector<uint8_t>visited(1u<<d);vector<uint32_t>todo;uint64_t covered=0;unsigned all=0,accepted=0;
 for(uint32_t s=0;s<(1u<<d);s++)if(!visited[s]){todo.clear();todo.push_back(s);visited[s]=1;bool even=true;
  for(size_t i=0;i<todo.size();i++){uint32_t x=todo[i];for(auto&r:a){uint32_t y=r.tab[0][x&255]|r.tab[1][x>>8&255]|r.tab[2][x>>16&255]|r.tab[3][x>>24&255];if(T)y^=r.flip;uint8_t color=1+((visited[x]-1)^(!T&&__builtin_parity(x&r.sign)));
   if(!visited[y]){visited[y]=color;todo.push_back(y);}else if(visited[y]!=color)even=false;}}
  covered+=todo.size();all++;if(T||even){accepted++;out<<s<<" "<<todo.size()<<"\n";}}
 if(covered!=(1u<<d))throw runtime_error("coverage");cerr<<"BFS all "<<all<<" accepted "<<accepted<<" covered "<<covered<<"\n";
}
