#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
using U=uint64_t;using T=std::vector<unsigned>;
U count(const T&r){unsigned n=r.size(),M=1u<<n;std::vector<U>d(M*n);for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;for(unsigned s=1;s<M;s++)for(unsigned v=0;v<n;v++)if(d[s*n+v]){unsigned a=r[v]&~s&(M-1);while(a){unsigned w=__builtin_ctz(a);a&=a-1;d[(s|(1u<<w))*n+w]+=d[s*n+v];}}U h=0;for(unsigned v=0;v<n;v++)h+=d[(M-1)*n+v];return h;}
T birth(T r,unsigned x){unsigned n=r.size();for(unsigned u=0;u<n;u++)if(u!=x)r[u]|=1u<<n;r.push_back(1u<<x);return r;}
U floorH(unsigned n){unsigned k=n/3;U a=1;for(unsigned i=0;i<k;i++)a*=5;if(n%3==0)return 3*a/5;if(n%3==1)return a;return 9*a/5;}
struct Hit{T rows;std::vector<unsigned> ports;std::vector<U> counts;};
struct Result{unsigned sources=0;std::map<unsigned,U> visits;std::map<U,Hit> hits;};
std::map<U,Result> results;
void dfs(U source,const T&root,const T&r,U h,std::vector<unsigned>&ports,std::vector<U>&counts){
 auto &o=results[source];o.visits[ports.size()]++;
 if(!ports.empty()&&(h==2*source-1||h==2*source+1)&&!o.hits.count(h))o.hits[h]={root,ports,counts};
 if(floorH(r.size()+1)>2*source+1)return;
 for(unsigned x=0;x<r.size();x++){T b=birth(r,x);U hh=count(b);if(hh<h)throw 1;if(hh>2*source+1)continue;ports.push_back(x);counts.push_back(hh);dfs(source,root,b,hh,ports,counts);counts.pop_back();ports.pop_back();}
}
int main(){unsigned n;U h;while(std::cin>>n>>h){T r(n);for(auto&v:r)std::cin>>v;long long a;for(unsigned i=0;i<n;i++)std::cin>>a;if(h!=31&&h!=43&&h!=61&&h!=73&&h!=113)continue;if(count(r)!=h)return 2;results[h].sources++;std::vector<unsigned>p;std::vector<U>c={h};dfs(h,r,r,h,p,c);}
 std::cout<<"{\"status\":\"PASS\",\"scope\":\"All single-port birth words from each strong source of the five prime counts, bounded by the target and inherited strong-order floor\",\"results\":[";bool sep=false;for(const auto &entry:results){if(sep)std::cout<<',';sep=true;const auto &o=entry.second;std::cout<<"{\"H\":"<<entry.first<<",\"sources\":"<<o.sources<<",\"accepted_nodes_by_depth\":{";bool s=false;for(auto v:o.visits){if(s)std::cout<<',';s=true;std::cout<<'"'<<v.first<<"\":"<<v.second;}std::cout<<"},\"hits\":[";s=false;for(auto a:o.hits){if(s)std::cout<<',';s=true;std::cout<<"{\"target\":"<<a.first<<",\"rows\":[";for(unsigned i=0;i<a.second.rows.size();i++){if(i)std::cout<<',';std::cout<<a.second.rows[i];}std::cout<<"],\"ports\":[";for(unsigned i=0;i<a.second.ports.size();i++){if(i)std::cout<<',';std::cout<<a.second.ports[i];}std::cout<<"],\"counts\":[";for(unsigned i=0;i<a.second.counts.size();i++){if(i)std::cout<<',';std::cout<<a.second.counts[i];}std::cout<<"]}";}std::cout<<"]}";}std::cout<<"]}\n";
}
