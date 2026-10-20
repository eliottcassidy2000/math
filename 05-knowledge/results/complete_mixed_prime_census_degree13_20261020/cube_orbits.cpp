#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
struct Action { uint32_t flip,sgn,weight,tab[4][256]{}; };
int main(int argc,char**argv){
 if(argc!=4)throw runtime_error("input output mode");
 ifstream in(argv[1]);ofstream out(argv[2]);int d,k;uint32_t N;in>>d>>k>>N;
 if(d>30||!in)throw runtime_error("invalid dimension");
 vector<Action> a(k);for(auto&r:a){in>>r.flip>>r.sgn>>r.weight;vector<unsigned>im(d);for(auto&v:im)in>>v;
  for(int b=0;b<4;b++)for(unsigned x=1;x<256;x++){unsigned z=__builtin_ctz(x),v=b*8+z;r.tab[b][x]=r.tab[b][x&(x-1)] | (v<(unsigned)d ? 1u<<im[v]:0);}}
 bool T=string(argv[3])=="T";vector<uint8_t>seen(1u<<d);unsigned total=0,accepted=0;
 for(uint32_t s=0;s<(1u<<d);s++)if(!seen[s]){
  bool even=true;uint32_t stab=0,size=0;
  for(auto&r:a){uint32_t y=r.tab[0][s&255]|r.tab[1][s>>8&255]|r.tab[2][s>>16&255]|r.tab[3][s>>24&255];if(T)y^=r.flip;
   if(!seen[y]){seen[y]=1;size++;}
   if(y==s){stab+=r.weight;if(!T&&__builtin_parity(s&r.sgn))even=false;}}
  if(size*uint64_t(stab)!=N)throw runtime_error("orbit stabilizer");total++;
  if(T||even){out<<s<<" "<<stab<<" "<<size<<"\n";accepted++;}}
 cerr<<"all "<<total<<" accepted "<<accepted<<" cube "<<(1u<<d)<<"\n";
}
