#include <array>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <vector>
using S=std::array<uint64_t,4>;
struct Hash{size_t operator()(S const& x)const{size_t h=0;for(auto a:x)h^=std::hash<uint64_t>{}(a)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
int main(){const uint64_t M=50000;std::unordered_set<S,Hash> cur{{0,0,0,1}};std::vector<bool> got(M+1);unsigned n=2;size_t total=0;while(!cur.empty()){
std::unordered_set<S,Hash> nxt;total+=cur.size();for(auto x:cur){auto[a,f,e,p]=x;auto h=a+2*e+p,r=a+f+2*e+2*p;got[h]=true;
auto offer=[&](S y){if(y[0]+2*y[2]+y[3]<=M)nxt.insert(y);};
offer({2*a+f+e+p,f,e+p,0});offer({f+p,0,a+e+p,e});
if(h+2*r<=M)offer({a+e,a+2*f+e+2*p,0,e+p});}
std::cout<<n<<" "<<cur.size()<<"\n";cur.swap(nxt);n++;}
std::cout<<"EXHAUSTED "<<n-1<<" TOTAL "<<total<<"\nMISSING";for(unsigned i=1;i<=M;i+=2)if(!got[i])std::cout<<" "<<i;std::cout<<"\n";}
