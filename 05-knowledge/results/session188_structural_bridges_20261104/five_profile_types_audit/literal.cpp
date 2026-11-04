#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
U count(std::vector<unsigned> r,unsigned a,unsigned c){unsigned n=r.size();r.push_back(1u<<(n+1));r.push_back(0);for(unsigned i=0;i<n;++i){if(a>>i&1)r[i]|=1u<<n;else r[n]|=1u<<i;if(c>>i&1)r[i]|=1u<<(n+1);else r[n+1]|=1u<<i;}n+=2;if(n>17)throw std::runtime_error("overflow bound");unsigned N=1u<<n;std::vector<U>d(size_t(N)*n);for(unsigned S=1;S<N;++S){unsigned bits=S;while(bits){unsigned v=__builtin_ctz(bits);bits&=bits-1;unsigned rest=S^(1u<<v);U z=rest?0:1;unsigned next=r[v]&rest;while(next){unsigned w=__builtin_ctz(next);next&=next-1;z+=d[size_t(rest)*n+w];}d[size_t(S)*n+v]=z;}}U z=0;for(unsigned i=0;i<n;++i)z+=d[size_t(N-1)*n+i];return z;}
int main(){unsigned n,k;while(std::cin>>n>>k){std::vector<unsigned>r(n);for(auto&v:r)std::cin>>v;while(k--){unsigned a,c;std::cin>>a>>c;std::cout<<count(r,a,c)<<'\n';}}}
