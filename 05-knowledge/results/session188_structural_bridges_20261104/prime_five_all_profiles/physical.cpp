#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
unsigned count(std::vector<unsigned> R,unsigned a,unsigned c){unsigned n=R.size();R.push_back(1u<<(n+1));R.push_back(0);for(unsigned i=0;i<n;++i){if(a>>i&1)R[i]|=1u<<n;else R[n]|=1u<<i;if(c>>i&1)R[i]|=1u<<(n+1);else R[n+1]|=1u<<i;}n+=2;if(n>22)throw std::runtime_error("physical cap22");unsigned N=1u<<n;std::vector<std::uint8_t>d(size_t(N)*n);for(unsigned S=1;S<N;++S){unsigned bits=S;while(bits){unsigned u=__builtin_ctz(bits);bits&=bits-1;unsigned rest=S^(1u<<u),z=rest?0:1,ps=rest&R[u];while(ps){unsigned v=__builtin_ctz(ps);ps&=ps-1;z+=d[size_t(rest)*n+v];}d[size_t(S)*n+u]=z&15;}}unsigned z=0;for(unsigned u=0;u<n;++u)z+=d[size_t(N-1)*n+u];return z&15;}
int main(){unsigned n,k;while(std::cin>>n>>k){std::vector<unsigned>T(n);for(auto&v:T)std::cin>>v;while(k--){unsigned a,c;std::cin>>a>>c;std::cout<<count(T,a,c)<<'\n';}}}
