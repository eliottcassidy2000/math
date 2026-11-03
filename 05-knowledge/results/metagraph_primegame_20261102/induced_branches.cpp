#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
using U=uint64_t;using T=std::vector<unsigned>;
std::vector<U> counts(const T&r){unsigned n=r.size(),M=1u<<n;std::vector<U>d(M*n),h(M);h[0]=1;for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;for(unsigned s=1;s<M;s++)for(unsigned v=0;v<n;v++)if(d[s*n+v]){h[s]+=d[s*n+v];unsigned a=r[v]&~s&(M-1);while(a){unsigned w=__builtin_ctz(a);a&=a-1;d[(s|(1u<<w))*n+w]+=d[s*n+v];}}return h;}
bool prime(unsigned p){for(unsigned d=2;d*d<=p;d++)if(p%d==0)return false;return p>=2;}
struct Hit{T rows;unsigned mask=0,added=99;};
int main(){unsigned n;U h,checks=0;std::map<std::pair<U,U>,Hit>hits;std::map<U,unsigned>sources;for(unsigned p=3;p<=61;p+=2)if(prime(p))for(U q:{2*p-1u,2*p+1u})hits[{p,q}]=Hit{};hits[{47,117}]=Hit{};
 while(std::cin>>n>>h){T r(n);for(auto&v:r)std::cin>>v;long long x;for(unsigned i=0;i<n;i++)std::cin>>x;auto hp=counts(r);if(hp.back()!=h)return 2;sources[h]++;
  for(unsigned mask=1;mask<hp.size()-1;mask++){
   checks++;auto key=std::make_pair(hp[mask],h);auto it=hits.find(key);if(it==hits.end())continue;unsigned added=n-__builtin_popcount(mask);if(added<it->second.added)it->second={r,mask,added};
  }
 }
 std::cout<<"{\"status\":\"PASS\",\"subset_checks\":"<<checks<<",\"branches\":[";bool sep=false;for(auto entry:hits){if(sep)std::cout<<',';sep=true;auto p=entry.first.first,q=entry.first.second;auto a=entry.second;std::cout<<"{\"p\":"<<p<<",\"q\":"<<q<<",\"target_strong_sources\":"<<sources[q]<<",\"embeds\":"<<(a.mask?"true":"false");if(a.mask){std::cout<<",\"added_vertices\":"<<a.added<<",\"source_subset\":"<<a.mask<<",\"target_rows\":[";for(unsigned i=0;i<a.rows.size();i++){if(i)std::cout<<',';std::cout<<a.rows[i];}std::cout<<']';}std::cout<<'}';}std::cout<<"]}\n";
}
