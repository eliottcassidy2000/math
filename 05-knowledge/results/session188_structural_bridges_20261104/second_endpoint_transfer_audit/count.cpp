#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned jobs;std::cin>>jobs;while(jobs--){unsigned n;std::cin>>n;if(n>18)throw std::runtime_error("order");std::vector<unsigned>R(n),P(n);for(auto&x:R)std::cin>>x;for(unsigned u=0;u<n;u++)for(unsigned v=0;v<n;v++)if(R[u]>>v&1)P[v]|=1u<<u;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n*3),s(std::size_t(N)*3);for(unsigned M=1;M<N;M++){unsigned vs=M;while(vs){unsigned b=vs&-vs;vs-=b;unsigned v=__builtin_ctz(b),rest=M^b;for(unsigned h=0;h<3;h++){U z=!rest&&h==0;if(rest){unsigned u=P[v]&rest;while(u){unsigned c=u&-u;u-=c;z+=d[(std::size_t(rest)*n+__builtin_ctz(c))*3+h];}if(h)z+=s[std::size_t(rest)*3+h-1];}d[(std::size_t(M)*n+v)*3+h]=z;s[std::size_t(M)*3+h]+=z;}}}for(unsigned h=0;h<3;h++)std::cout<<s[std::size_t(N-1)*3+h]<<' ';std::cout<<'\n';}}
