#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
using U=std::uint64_t;
std::pair<U,U> count(const std::vector<unsigned>&old,unsigned a,unsigned c){
 unsigned n=old.size(),m=n+2;if(m>17)throw std::runtime_error("literal order exceeds17");std::vector<unsigned>R=old;R.push_back(1u<<(n+1));R.push_back(0);
 for(unsigned u=0;u<n;++u){if(a>>u&1u)R[u]|=1u<<n;else R[n]|=1u<<u;if(c>>u&1u)R[u]|=1u<<(n+1);else R[n+1]|=1u<<u;}
 static std::map<std::vector<unsigned>,std::pair<U,U>> cache; auto found=cache.find(R); if(found!=cache.end()) return found->second;
 unsigned N=1u<<m;std::vector<U>d(std::size_t(N)*m);for(unsigned u=0;u<m;++u)d[std::size_t(1u<<u)*m+u]=1;
 for(unsigned S=1;S<N;++S){unsigned left=S;while(left){unsigned u=__builtin_ctz(left);left&=left-1;U z=d[std::size_t(S)*m+u];if(!z)continue;unsigned next=R[u]&~S;while(next){unsigned b=next&-next;next-=b;unsigned v=__builtin_ctz(b);d[std::size_t(S|b)*m+v]+=z;}}}
 std::vector<U> H(N); for(unsigned S=1;S<N;++S)for(unsigned u=0;u<m;++u)H[S]+=d[std::size_t(S)*m+u]; U q2=0;for(unsigned S=1;S<N-1;++S)q2+=H[S]*H[(N-1)^S]; auto answer=std::make_pair(H[N-1],q2);cache.emplace(R,answer);return answer;
}
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n,k;std::cin>>n>>k;std::vector<unsigned>T(n),V(n);for(auto&x:T)std::cin>>x;for(auto&x:V)std::cin>>x;while(k--){unsigned a,c;std::cin>>a>>c;auto A=count(T,a,c),B=count(V,a,c); std::cout<<A.first<<' '<<B.first<<' '<<A.second<<' '<<B.second<<'\n';}}}
