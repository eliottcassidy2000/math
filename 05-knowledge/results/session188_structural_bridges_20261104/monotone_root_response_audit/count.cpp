#include <cstdint>
#include <vector>
#include <iostream>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned count;std::cin>>count;while(count--){unsigned n;std::cin>>n;if(n>18)throw std::runtime_error("physical order exceeds18");std::vector<unsigned>a(n);for(auto&x:a)std::cin>>x;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n);for(unsigned S=1;S<N;++S){unsigned rem=S;while(rem){unsigned u=__builtin_ctz(rem);rem&=rem-1;unsigned prev=S^(1u<<u);U z=0;if(!prev)z=1;else{unsigned next=a[u]&prev;while(next){unsigned v=__builtin_ctz(next);next&=next-1;z+=d[std::size_t(prev)*n+v];}}d[std::size_t(S)*n+u]=z;}}U answer=0;for(unsigned u=0;u<n;++u)answer+=d[std::size_t(N-1)*n+u];std::cout<<answer<<'\n';}}
