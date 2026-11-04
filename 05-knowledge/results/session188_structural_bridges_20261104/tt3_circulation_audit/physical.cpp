#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned n;while(std::cin>>n){if(n>20)throw std::runtime_error("literal exact cap20");std::vector<unsigned>R(n);for(auto&v:R)std::cin>>v;unsigned N=1u<<n;std::vector<U>d(size_t(N)*n);for(unsigned M=1;M<N;++M){unsigned bits=M;while(bits){unsigned u=__builtin_ctz(bits);bits&=bits-1;unsigned S=M^(1u<<u),vs=R[u]&S;U z=S?0:1;while(vs){unsigned v=__builtin_ctz(vs);vs&=vs-1;z+=d[size_t(S)*n+v];}d[size_t(M)*n+u]=z;}}U z=0;for(unsigned u=0;u<n;++u)z+=d[size_t(N-1)*n+u];std::cout<<z<<'\n';}}
