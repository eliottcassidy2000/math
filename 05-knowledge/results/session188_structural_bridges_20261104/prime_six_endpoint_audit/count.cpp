#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
using U=std::uint64_t;
std::vector<U> count(const std::vector<unsigned>&out){unsigned n=out.size();if(n>17)throw std::runtime_error("order");unsigned N=1u<<n;std::vector<unsigned>in(n);for(unsigned i=0;i<n;++i)for(unsigned j=0;j<n;++j)if(out[i]>>j&1)in[j]|=1u<<i;std::vector<U>d(std::size_t(N)*n*3),total(std::size_t(N)*3);for(unsigned M=1;M<N;++M){unsigned bits=M;while(bits){unsigned v=__builtin_ctz(bits);bits&=bits-1;unsigned R=M^(1u<<v);for(unsigned h=0;h<3;++h){U z=R?(h?total[std::size_t(R)*3+h-1]:0):(h==0);unsigned pre=R&in[v];while(pre){unsigned u=__builtin_ctz(pre);pre&=pre-1;z+=d[(std::size_t(R)*n+u)*3+h];}d[(std::size_t(M)*n+v)*3+h]=z;total[std::size_t(M)*3+h]+=z;}}}return {total[std::size_t(N-1)*3],total[std::size_t(N-1)*3+1],total[std::size_t(N-1)*3+2]};}
int main(){unsigned cases;std::cin>>cases;std::map<std::vector<unsigned>,std::vector<U>>cache;while(cases--){unsigned n;std::cin>>n;std::vector<unsigned>R(n);for(auto&r:R)std::cin>>r;auto it=cache.find(R);if(it==cache.end())it=cache.emplace(R,count(R)).first;for(auto z:it->second)std::cout<<z<<' ';std::cout<<'\n';}}
