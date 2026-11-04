#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned cases;if(!(std::cin>>cases))return 1;while(cases--){unsigned n;std::cin>>n;if(n>18)throw std::runtime_error("order");std::vector<unsigned>out(n),in(n);for(auto &r:out)std::cin>>r;for(unsigned u=0;u<n;++u)for(unsigned v=0;v<n;++v)if(out[u]>>v&1)in[v]|=1u<<u;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n);for(unsigned M=1;M<N;++M){unsigned z=M;while(z){unsigned v=__builtin_ctz(z);z&=z-1;unsigned rest=M^(1u<<v);U value=rest?0:1;unsigned p=rest&in[v];while(p){unsigned u=__builtin_ctz(p);p&=p-1;value+=d[std::size_t(rest)*n+u];}d[std::size_t(M)*n+v]=value;}}U ans=0;for(unsigned u=0;u<n;++u)ans+=d[std::size_t(N-1)*n+u];std::cout<<ans<<'\n';}}
