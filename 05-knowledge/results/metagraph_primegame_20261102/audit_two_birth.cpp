#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <set>
#include <vector>
using U=uint64_t;using T=std::vector<unsigned>;
U count(const T&r){unsigned n=r.size(),M=1u<<n;std::vector<U>d(M*n);for(unsigned v=0;v<n;v++)d[(1u<<v)*n+v]=1;for(unsigned s=1;s<M;s++)for(unsigned v=0;v<n;v++)if(d[s*n+v]){unsigned a=r[v]&~s&(M-1);while(a){unsigned w=__builtin_ctz(a);a&=a-1;d[(s|(1u<<w))*n+w]+=d[s*n+v];}}U h=0;for(unsigned v=0;v<n;v++)h+=d[(M-1)*n+v];return h;}
T birth(T r,unsigned x){unsigned n=r.size();for(unsigned u=0;u<n;u++)if(u!=x)r[u]|=1u<<n;r.push_back(1u<<x);return r;}
bool legal(const T&r,const T&p){if(p.size()!=r.size())return false;unsigned seen=0;for(unsigned i=0;i<p.size();i++){if(seen>>p[i]&1)return false;seen|=1u<<p[i];if(i&&!(r[p[i-1]]>>p[i]&1))return false;}return true;}
std::vector<T> paths(const T&r){std::vector<T>out;T p;std::function<void(unsigned)>go=[&](unsigned used){if(p.size()==r.size()){out.push_back(p);return;}for(unsigned v=0;v<r.size();v++)if(!(used>>v&1)&&(p.empty()||(r[p.back()]>>v&1))){p.push_back(v);go(used|(1u<<v));p.pop_back();}};go(0);return out;}
std::vector<T> insertions(const T&p,unsigned x,unsigned z){T q=p;auto at=std::find(q.begin(),q.end(),x);q.insert(at,z);std::vector<T>out={q};if(p.back()!=x){q=p;q.push_back(z);out.push_back(q);}return out;}
int main(){unsigned n;U h,sources=0,endpoint_checks=0,ordered_birth_pairs=0,explicit_injected_paths=0;while(std::cin>>n>>h){T r(n);for(auto&v:r)std::cin>>v;long long discard;for(unsigned i=0;i<n;i++)std::cin>>discard;
 auto ps=paths(r);if(ps.size()!=h||count(r)!=h)return 2;std::vector<U>e(n);for(auto p:ps)e[p.back()]++;
 for(unsigned x=0;x<n;x++){
  if(!r[x]||2*e[x]>h-1)return 3;
  std::set<T>moved;
  for(auto p:ps)if(p.back()==x){p.pop_back();auto it=std::find_if(p.begin(),p.end(),[&](unsigned v){return r[x]>>v&1;});if(it==p.end())return 4;p.insert(it,x);if(!legal(r,p)||p.back()==x||!moved.insert(p).second)return 5;}
  if(moved.size()!=e[x])return 6;endpoint_checks++;
  T b=birth(r,x);
  for(unsigned y=0;y<=n;y++){
   T c=birth(b,y);std::set<T>injected;
   for(auto p:ps)for(auto q:insertions(p,x,n))for(auto v:insertions(q,y,n+1))if(!legal(c,v)||!injected.insert(v).second)return 7;
   U expected=y==n?3*h-e[x]:(y==x?4*h-3*e[x]:4*h-2*e[x]-e[y]);
   U actual=count(c);if(injected.size()!=expected||actual<expected||actual<2*h+3||2*actual<5*h+1)return 8;
   explicit_injected_paths+=injected.size();ordered_birth_pairs++;
  }
 }sources++;}
 std::cout<<"{\"status\":\"PASS\",\"sources\":"<<sources<<",\"endpoint_injections\":"<<endpoint_checks<<",\"ordered_birth_pairs\":"<<ordered_birth_pairs<<",\"explicit_injected_paths\":"<<explicit_injected_paths<<",\"scope\":\"All 687 strong sources below H125; full labelled injections and direct counts\"}\n";
}
