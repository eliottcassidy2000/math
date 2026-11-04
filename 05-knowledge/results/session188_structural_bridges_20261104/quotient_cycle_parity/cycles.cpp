#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
using U=std::uint64_t;
U cycles(const std::vector<unsigned>&R){
 unsigned n=R.size();if(n>18 || n<3)throw std::runtime_error("cycle order outside3..18");
 static std::map<std::vector<unsigned>,U> cache;auto it=cache.find(R);if(it!=cache.end())return it->second;
 unsigned root=n-2,N=1u<<n,b=1u<<root;std::vector<U>d(std::size_t(N)*n);d[std::size_t(b)*n+root]=1;
 for(unsigned S=1;S<N;++S){if(!(S&b))continue;unsigned rem=S;while(rem){unsigned u=__builtin_ctz(rem);rem&=rem-1;U z=d[std::size_t(S)*n+u];if(!z)continue;unsigned next=R[u]&~S;while(next){unsigned t=next&-next;next-=t;unsigned v=__builtin_ctz(t);d[std::size_t(S|t)*n+v]+=z;}}}
 U out=0;for(unsigned u=0;u<n;++u)if(R[u]&b)out+=d[std::size_t(N-1)*n+u];cache[R]=out;return out;
}
int main(){unsigned k;std::cin>>k;while(k--){unsigned n;std::cin>>n;std::vector<unsigned>R(n);for(auto&x:R)std::cin>>x;std::cout<<cycles(R)<<'\n';}}
