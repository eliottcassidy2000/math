#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n,k,a,c;std::cin>>n>>k;std::vector<unsigned>R(n);for(auto&x:R)std::cin>>x;std::cin>>a>>c;unsigned m=n+2;if(m>14||k>6)throw std::runtime_error("audit bound");R.resize(m);R[n]=1u<<(n+1);for(unsigned u=0;u<n;++u){if(a>>u&1)R[u]|=1u<<n;else R[n]|=1u<<u;if(c>>u&1)R[u]|=1u<<(n+1);else R[n+1]|=1u<<u;}std::vector<unsigned>pred(m);for(unsigned u=0;u<m;++u)for(unsigned v=0;v<m;++v)if(R[u]>>v&1)pred[v]|=1u<<u;unsigned N=1u<<m;std::vector<U>D(size_t(N)*m*(k+1));auto at=[&](unsigned M,unsigned v,unsigned h)->U&{return D[(size_t(M)*m+v)*(k+1)+h];};for(unsigned M=1;M<N;++M){unsigned vertices=M;while(vertices){unsigned v=__builtin_ctz(vertices);vertices&=vertices-1;unsigned S=M^(1u<<v);if(!S){at(M,v,v<k)=1;continue;}unsigned p=pred[v]&S;while(p){unsigned u=__builtin_ctz(p);p&=p-1;unsigned delta=(v<k&&u>=k);for(unsigned h=delta;h<=k;++h)at(M,v,h)+=at(S,u,h-delta);}}}for(unsigned h=0;h<=k;++h){U z=0;for(unsigned v=0;v<m;++v)z+=at(N-1,v,h);if(h)std::cout<<' ';std::cout<<z;}std::cout<<'\n';}}
