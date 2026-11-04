#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n,K;std::cin>>n>>K;if(n>18||K<1||K>5||(K>1&&n>16))throw std::runtime_error("unsupported physical bound");std::vector<unsigned>A(n);for(auto&x:A)std::cin>>x;unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n),H(N);for(unsigned S=1;S<N;++S){unsigned rem=S;while(rem){unsigned v=__builtin_ctz(rem);rem&=rem-1;unsigned prev=S^(1u<<v);U z=prev?0:1;unsigned next=A[v]&prev;while(next){unsigned w=__builtin_ctz(next);next&=next-1;z+=d[std::size_t(prev)*n+w];}d[std::size_t(S)*n+v]=z;H[S]+=z;}}std::cout<<H[N-1];if(K>1){std::vector<U>part(std::size_t(K+1)*N);part[0]=1;for(unsigned S=1;S<N;++S){unsigned low=S&-S,rest=S^low,sub=rest;while(true){unsigned block=sub|low,tail=S^block;for(unsigned k=1;k<=K;++k)part[std::size_t(k)*N+S]+=H[block]*part[std::size_t(k-1)*N+tail];if(!sub)break;sub=(sub-1)&rest;}}U fac=1;for(unsigned k=1;k<=K;++k){fac*=k;if(k>1)std::cout<<' '<<fac*part[std::size_t(k)*N+N-1];}}std::cout<<'\n';}return std::cin.fail()?1:0;}
