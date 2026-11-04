#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
int main(){unsigned jobs;std::cin>>jobs;while(jobs--){unsigned n;std::cin>>n;if(n>18)throw std::runtime_error("order");std::vector<unsigned>R(n),P(n);for(auto&v:R)std::cin>>v;for(unsigned u=0;u<n;u++)for(unsigned v=0;v<n;v++)if((R[u]>>v)&1)P[v]|=1u<<u;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n),h(N);for(unsigned S=1;S<N;S++){unsigned bs=S;while(bs){unsigned b=bs&-bs;bs-=b;unsigned v=__builtin_ctz(b),rest=S^b;U z=0;if(!rest)z=1;else{unsigned pre=P[v]&rest;while(pre){unsigned q=pre&-pre;pre-=q;z+=d[std::size_t(rest)*n+__builtin_ctz(q)];}}d[std::size_t(S)*n+v]=z;h[S]+=z;}}U two=0;for(unsigned S=1;S<N-1;S++)two+=h[S]*h[(N-1)^S];std::cout<<h[N-1]<<' '<<two<<'\n';}}
