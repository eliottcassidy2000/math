#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using U=std::uint64_t;
U physical(const std::vector<unsigned>& a){
 const unsigned n=a.size();if(n>17)throw std::runtime_error("order exceeds17");
 const unsigned N=1u<<n;std::vector<U>d(std::size_t(N)*n);
 for(unsigned S=1;S<N;++S){unsigned vertices=S;while(vertices){unsigned v=__builtin_ctz(vertices);vertices&=vertices-1;unsigned prev=S^(1u<<v);U z=0;
  if(!prev)z=1;else{unsigned destinations=a[v]&prev;while(destinations){unsigned w=__builtin_ctz(destinations);destinations&=destinations-1;z+=d[std::size_t(prev)*n+w];}}
  d[std::size_t(S)*n+v]=z;
 }}
 U out=0;for(unsigned v=0;v<n;++v)out+=d[std::size_t(N-1)*n+v];return out;
}
U run(const std::array<unsigned,5>&Q,unsigned a,unsigned c,const std::array<unsigned,5>&h){
 std::array<unsigned,7>R{};for(unsigned u=0;u<5;++u)R[u]=Q[u];R[5]=1u<<6;
 for(unsigned u=0;u<5;++u){if(a>>u&1u)R[u]|=1u<<5;else R[5]|=1u<<u;if(c>>u&1u)R[u]|=1u<<6;else R[6]|=1u<<u;}
 std::vector<U>memo(4096*7,std::numeric_limits<U>::max());
 std::function<U(unsigned,unsigned)>f=[&](unsigned S,unsigned last)->U{
  U&out=memo[S*7+last];if(out!=std::numeric_limits<U>::max())return out;
  bool terminal=(S>>10)==3;for(unsigned i=0;i<5;++i)terminal=terminal&&((S>>(2*i))&3u);
  out=terminal?1:0;unsigned next=R[last];
  while(next){unsigned v=__builtin_ctz(next);next&=next-1;unsigned T,weight=1;
   if(v<5){unsigned k=(S>>(2*v))&3u;if(k==3)continue;T=S+(1u<<(2*v));weight=k==0?h[v]:(k==1?6/h[v]:1);}
   else{unsigned bit=1u<<(10+v-5);if(S&bit)continue;T=S|bit;}
   out+=weight*f(T,v);
  }return out;
 };
 U out=f(1u<<10,5)+f(1u<<11,6);for(unsigned v=0;v<5;++v)out+=h[v]*f(1u<<(2*v),v);return out;
}
int main(){unsigned jobs;if(!(std::cin>>jobs))return 1;for(unsigned t=0;t<jobs;++t){char kind;std::cin>>kind;if(kind=='P'){unsigned n;std::cin>>n;std::vector<unsigned>a(n);for(auto&v:a)std::cin>>v;std::cout<<physical(a)<<'\n';}else if(kind=='R'){std::array<unsigned,5>Q,h;unsigned a,c;for(auto&v:Q)std::cin>>v;std::cin>>a>>c;for(auto&v:h)std::cin>>v;std::cout<<run(Q,a,c,h)<<'\n';}else throw std::runtime_error("invalid mode");}}
