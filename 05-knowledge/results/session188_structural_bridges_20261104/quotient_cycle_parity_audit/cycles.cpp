#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned jobs;std::cin>>jobs;while(jobs--){unsigned n;std::cin>>n;if(n<3||n>18)throw std::runtime_error("order");std::vector<unsigned>R(n),P(n);for(auto&x:R)std::cin>>x;for(unsigned u=0;u<n;u++)for(unsigned v=0;v<n;v++)if(R[u]>>v&1)P[v]|=1u<<u;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n);d[n]=1;for(unsigned M=3;M<N;M+=2){unsigned vs=M&~1u;while(vs){unsigned b=vs&-vs;vs-=b;unsigned v=__builtin_ctz(b),rest=M^b,pre=P[v]&rest;U z=0;while(pre){unsigned t=pre&-pre;pre-=t;z+=d[std::size_t(rest)*n+__builtin_ctz(t)];}d[std::size_t(M)*n+v]=z;}}U z=0;for(unsigned v=1;v<n;v++)if(R[v]&1)z+=d[std::size_t(N-1)*n+v];std::cout<<z<<'\n';}}
