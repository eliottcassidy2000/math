#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
using S=std::array<uint64_t,4>;
struct Hash {size_t operator()(S const&x)const{size_t h=0;for(auto a:x)h^=std::hash<uint64_t>{}(a)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
int main(){constexpr uint64_t M=50000; std::unordered_map<S,uint64_t,Hash> cur{{{1,2,0,1},0}};std::vector<uint64_t> words(M+1,UINT64_MAX);unsigned n=2;uint64_t total=0;
while(!cur.empty()){std::unordered_map<S,uint64_t,Hash> nxt;total+=cur.size();for(auto const&item:cur){auto[h,r,e,d]=item.first;auto w=item.second; if(words[h]==UINT64_MAX||w<words[h])words[h]=w;
auto offer=[&](S y,unsigned letter){if(n>=32){std::cerr<<"Word encoding limit reached\n";std::exit(2);}auto nw=4*w+letter;auto[it,fresh]=nxt.emplace(y,nw);if(!fresh&&nw<it->second)it->second=nw;};
if(h+r-e<=M){offer({h+r-e,2*r-d,d,d},2);offer({h+r-e,h+r,h-e,h},3);}if(h+2*r<=M)offer({h,2*r,0,d},1);}
std::cout<<n<<" "<<cur.size()<<"\n";cur.swap(nxt);++n;}
std::ofstream f("WITNESSES.tsv");unsigned count=0;for(unsigned h=1;h<=M;h+=2)if(words[h]!=UINT64_MAX){auto w=words[h];std::string s;while(w){s.push_back("?OZJ"[w%4]);w/=4;}std::reverse(s.begin(),s.end());f<<h<<"\t"<<(s.empty()?"-":s)<<"\n";++count;}
std::cout<<"EXHAUSTED "<<n-1<<" TOTAL "<<total<<" WITNESSES "<<count<<"\nMISSING";for(unsigned h=1;h<=M;h+=2)if(words[h]==UINT64_MAX)std::cout<<" "<<h;std::cout<<"\n";
}
