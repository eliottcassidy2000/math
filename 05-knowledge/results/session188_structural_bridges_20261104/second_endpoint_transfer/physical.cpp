#include <cstdint>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
using U=std::uint64_t;
std::vector<U> count(const std::vector<unsigned>&R){unsigned n=R.size();if(n>18)throw std::runtime_error("physical order");unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n*3);auto at=[&](unsigned M,unsigned v,unsigned h)->U&{return d[(std::size_t(M)*n+v)*3+h];};for(unsigned v=0;v<n;++v)at(1u<<v,v,0)=1;for(unsigned M=1;M<N;++M){unsigned old=M;while(old){unsigned v=__builtin_ctz(old);old&=old-1;for(unsigned h=0;h<3;++h){U z=at(M,v,h);if(!z)continue;unsigned next=(N-1)^M;while(next){unsigned w=__builtin_ctz(next);next&=next-1;unsigned MM=M|(1u<<w);if(R[v]>>w&1)at(MM,w,h)+=z;if(h<2)at(MM,w,h+1)+=z;}}}}std::vector<U>q(3);for(unsigned v=0;v<n;++v)for(unsigned h=0;h<3;++h)q[h]+=at(N-1,v,h);return q;}
int main(){unsigned cases;std::cin>>cases;std::map<std::vector<unsigned>,std::vector<U>>cache;while(cases--){unsigned n;std::cin>>n;std::vector<unsigned>R(n);for(auto&r:R)std::cin>>r;auto it=cache.find(R);if(it==cache.end())it=cache.emplace(R,count(R)).first;for(auto q:it->second)std::cout<<q<<' ';std::cout<<'\n';}}
