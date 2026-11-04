#include <cstdint>
#include <iostream>
#include <vector>
int main(){unsigned n;while(std::cin>>n){if(n>20)return 3;std::vector<unsigned>row(n);for(auto &r:row)std::cin>>r;const unsigned N=1u<<n;std::vector<uint64_t>d(uint64_t(N)*n);for(unsigned v=0;v<n;v++)d[uint64_t(1u<<v)*n+v]=1;for(unsigned mask=1;mask<N;mask++)for(unsigned ends=mask;ends;ends&=ends-1){unsigned v=__builtin_ctz(ends);uint64_t ways=d[uint64_t(mask)*n+v];for(unsigned b=row[v]&~mask;b;b&=b-1){unsigned u=__builtin_ctz(b);d[uint64_t(mask|(1u<<u))*n+u]+=ways;}}uint64_t answer=0;for(unsigned v=0;v<n;v++)answer+=d[uint64_t(N-1)*n+v];std::cout<<answer<<'\n';}}
