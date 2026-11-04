#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <map>
using U=std::uint64_t;
std::vector<U> count(const std::vector<unsigned>&R,bool covers){
 unsigned n=R.size(),N=1u<<n;if(n>20||(covers&&n>12))throw std::runtime_error("bound");std::vector<U>d(std::size_t(N)*n),H(N);
 for(unsigned u=0;u<n;++u)d[std::size_t(1u<<u)*n+u]=1;
 for(unsigned M=1;M<N;++M){unsigned bits=M;while(bits){unsigned u=__builtin_ctz(bits);bits&=bits-1;U z=d[std::size_t(M)*n+u];H[M]+=z;if(!z)continue;unsigned next=R[u]&~M;while(next){unsigned v=__builtin_ctz(next);next&=next-1;d[std::size_t(M|(1u<<v))*n+v]+=z;}}}
 if(!covers){return{H.back()};}std::vector<U>B(N);for(unsigned M=1;M<N;++M)for(unsigned S=(M-1)&M;S;S=(S-1)&M)B[M]+=H[S]*H[M^S];U q3=0;for(unsigned S=1;S<N-1;++S)q3+=B[S]*H[(N-1)^S];return{H.back(),B.back(),q3};
}
int main(){unsigned cases;std::cin>>cases;std::map<std::vector<unsigned>,std::vector<U>>cache;while(cases--){unsigned n,cover;std::cin>>n>>cover;std::vector<unsigned>R(n);for(auto&r:R)std::cin>>r;auto key=R;key.push_back(cover);auto it=cache.find(key);if(it==cache.end())it=cache.emplace(key,count(R,cover)).first;for(auto v:it->second)std::cout<<v<<' ';std::cout<<'\n';}}
