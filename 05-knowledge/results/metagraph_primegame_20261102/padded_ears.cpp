#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using U=uint64_t;using T=std::vector<unsigned>;
U count(const T&r){unsigned n=r.size(),M=1u<<n;std::vector<U>d(M*n);for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;for(unsigned s=1;s<M;s++)for(unsigned v=0;v<n;v++)if(d[s*n+v]){unsigned a=r[v]&~s&(M-1);while(a){unsigned w=__builtin_ctz(a);a&=a-1;d[(s|(1u<<w))*n+w]+=d[s*n+v];}}U h=0;for(unsigned v=0;v<n;v++)h+=d[(M-1)*n+v];return h;}
T padded(const T&r,unsigned a,unsigned b){unsigned n=r.size(),N=a+n+b;T out(N);for(unsigned u=0;u<N;u++)for(unsigned v=u+1;v<N;v++){if(u>=a&&v<a+n){if(r[u-a]>>(v-a)&1)out[u]|=1u<<v;else out[v]|=1u<<u;}else out[u]|=1u<<v;}return out;}
T ear(T r,unsigned mask){unsigned n=r.size();for(unsigned u=0;u<n;u++)if(!(mask>>u&1))r[u]|=1u<<n;r.push_back(mask);return r;}
struct Hit{U h,target;unsigned a,b,mask,ports;T core;};
int main(){unsigned n;U h,checks=0,sources=0;std::set<U>dead={31,43,55,61,63,73,85,87,113};std::map<std::pair<U,int>,Hit>hits;while(std::cin>>n>>h){T r(n);for(auto&v:r)std::cin>>v;long long x;for(unsigned i=0;i<n;i++)std::cin>>x;if(!dead.count(h))continue;
 for(unsigned a=0;a<=2;a++)for(unsigned b=0;b+a<=2;b++){
  T p=padded(r,a,b);if(count(p)!=h)return 2;sources++;
  for(unsigned mask=1;mask<(1u<<p.size())-1;mask++){
   U target=count(ear(p,mask));checks++;if(target!=2*h-1&&target!=2*h+1)continue;
   unsigned ports=__builtin_popcount(mask);if(ports==1)return 3;int sign=target==2*h-1?-1:1;auto key=std::make_pair(h,sign);
   if(!hits.count(key)||ports<hits[key].ports)hits[key]={h,target,a,b,mask,ports,r};
  }
 }
}
 std::cout<<"{\"status\":\"PASS\",\"scope\":\"The nine one-port-dead counts; all strong cores below125, all prefix/suffix lengths with sum at most2, all fresh-vertex cuts\",\"sources\":"<<sources<<",\"checks\":"<<checks<<",\"hits\":[";bool sep=false;for(auto pair:hits){auto a=pair.second;if(sep)std::cout<<',';sep=true;std::cout<<"{\"H\":"<<a.h<<",\"target\":"<<a.target<<",\"prefix\":"<<a.a<<",\"suffix\":"<<a.b<<",\"mask\":"<<a.mask<<",\"ports\":"<<a.ports<<",\"core\":[";for(unsigned i=0;i<a.core.size();i++){if(i)std::cout<<',';std::cout<<a.core[i];}std::cout<<"]}";}std::cout<<"]}\n";
}
