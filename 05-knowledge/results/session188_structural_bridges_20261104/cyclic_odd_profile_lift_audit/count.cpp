#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using U=std::uint64_t;
int main(){unsigned cases;std::cin>>cases;while(cases--){unsigned n,a,c;std::cin>>n;std::vector<unsigned>R(n);for(auto&x:R)std::cin>>x;std::cin>>a>>c;unsigned m=n+2;if(m>17)throw std::runtime_error("order bound");R.resize(m);R[n]=1u<<(n+1);for(unsigned u=0;u<n;++u){if(a>>u&1)R[u]|=1u<<n;else R[n]|=1u<<u;if(c>>u&1)R[u]|=1u<<(n+1);else R[n+1]|=1u<<u;}std::vector<unsigned>pred(m);for(unsigned u=0;u<m;++u)for(unsigned v=0;v<m;++v)if(R[u]>>v&1)pred[v]|=1u<<u;unsigned N=1u<<m;std::vector<U>D(size_t(N)*m),H(N);for(unsigned M=1;M<N;++M){unsigned vertices=M;while(vertices){unsigned v=__builtin_ctz(vertices);vertices&=vertices-1;unsigned S=M^(1u<<v);U z=0;if(!S)z=1;else{unsigned p=pred[v]&S;while(p){unsigned w=__builtin_ctz(p);p&=p-1;z+=D[size_t(S)*m+w];}}D[size_t(M)*m+v]=z;H[M]+=z;}}U pair=0;for(unsigned M=1;M<N-1;++M)if(M&1)pair+=2*H[M]*H[(N-1)^M];std::cout<<H[N-1]<<' '<<pair<<'\n';}}
