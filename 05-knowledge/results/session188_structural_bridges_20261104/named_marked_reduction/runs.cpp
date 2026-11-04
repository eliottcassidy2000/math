#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
std::vector<U> counts(std::vector<unsigned> rows,unsigned a,unsigned c,unsigned m){unsigned n=rows.size();rows.push_back(1u<<(n+1));rows.push_back(0);for(unsigned i=0;i<n;++i){if(a>>i&1)rows[i]|=1u<<n;else rows[n]|=1u<<i;if(c>>i&1)rows[i]|=1u<<(n+1);else rows[n+1]|=1u<<i;}n+=2;if(n>14)throw std::runtime_error("overflow cap");unsigned N=1u<<n,K=m+1;std::vector<U>d(size_t(N)*n*K);auto at=[&](unsigned S,unsigned v,unsigned k)->U&{return d[(size_t(S)*n+v)*K+k];};for(unsigned i=0;i<n;++i)at(1u<<i,i,i<m)=1;for(unsigned S=1;S<N;++S){unsigned last=S;while(last){unsigned u=__builtin_ctz(last);last&=last-1;unsigned next=rows[u]&~S;while(next){unsigned v=__builtin_ctz(next);next&=next-1;unsigned add=(v<m&&u>=m);for(unsigned k=0;k+add<K;++k)at(S|(1u<<v),v,k+add)+=at(S,u,k);}}}std::vector<U>out(K);for(unsigned u=0;u<n;++u)for(unsigned k=0;k<K;++k)out[k]+=at(N-1,u,k);return out;}
int main(){unsigned n,m,k;while(std::cin>>n>>m>>k){std::vector<unsigned>r(n);for(auto&v:r)std::cin>>v;while(k--){unsigned a,c;std::cin>>a>>c;auto vv=counts(r,a,c,m);for(auto v:vv)std::cout<<v<<' ';std::cout<<'\n';}}}
