#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned n;while(std::cin>>n){if(n>14)throw std::runtime_error("exact cover cap14");std::vector<unsigned>R(n),pred(n);for(auto&v:R)std::cin>>v;for(unsigned u=0;u<n;++u)for(unsigned v=0;v<n;++v)if(R[u]>>v&1)pred[v]|=1u<<u;unsigned N=1u<<n;std::vector<U>d(size_t(N)*n*3),tot(size_t(N)*3);for(unsigned M=1;M<N;++M){unsigned vs=M;while(vs){unsigned v=__builtin_ctz(vs);vs&=vs-1;unsigned S=M^(1u<<v);U z[3]={0,0,0};if(!S)z[0]=1;else{z[1]=tot[size_t(S)*3];z[2]=tot[size_t(S)*3+1];unsigned ps=pred[v]&S;while(ps){unsigned u=__builtin_ctz(ps);ps&=ps-1;for(unsigned k=0;k<3;++k)z[k]+=d[(size_t(S)*n+u)*3+k];}}for(unsigned k=0;k<3;++k){d[(size_t(M)*n+v)*3+k]=z[k];tot[size_t(M)*3+k]+=z[k];}}}for(unsigned k=0;k<3;++k)std::cout<<tot[size_t(N-1)*3+k]<<' ';std::cout<<'\n';}}
