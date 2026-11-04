#include <algorithm>
#include <array>
#include <iostream>
#include <cstdint>
int main(){int nq;std::cin>>nq;while(nq--){std::array<unsigned,8>Q;for(auto&r:Q)std::cin>>r;bool prime=true;for(unsigned m=1;m<255;++m)if(__builtin_popcount(m)>1){bool module=true;for(int i=0;i<8;++i)if(!(m>>i&1)){auto x=Q[i]&m;if(x!=0&&x!=m){module=false;break;}}if(module){prime=false;break;}}std::array<int,8>p={0,1,2,3,4,5,6,7};std::array<uint64_t,8>S{},E{};do{bool good=true;for(int i=0;i<7;++i)if(!(Q[p[i]]>>p[i+1]&1)){good=false;break;}if(good){S[p[0]]++;E[p[7]]++;}}while(std::next_permutation(p.begin(),p.end()));std::cout<<prime;for(auto x:S)std::cout<<' '<<x;for(auto x:E)std::cout<<' '<<x;std::cout<<'\n';}}
