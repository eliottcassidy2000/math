// Independent deletion recurrence. Unlike search.cpp, this never propagates
// a partial path forward: each cell is reconstructed from its smaller subset.
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace std;
int64_t add(int64_t a,int64_t b){if(b>numeric_limits<int64_t>::max()-a)throw overflow_error("sum overflow");return a+b;}
int64_t mul(int64_t a,int64_t b){if(a && b>numeric_limits<int64_t>::max()/a)throw overflow_error("product overflow");return a*b;}
int main(){int n;while(cin>>n){if(n<1||n>21)return 2;unsigned full=(1u<<n)-1;
 vector<unsigned>t(n),in(n);for(auto &v:t)cin>>v;
 for(int v=0;v<n;v++){if(t[v]&(1u<<v)||t[v]&~full)return 3;
  for(int u=0;u<v;u++)if(((t[v]>>u)&1)==((t[u]>>v)&1))return 4;
  in[v]=full^t[v]^(1u<<v);}
 vector<int64_t>d((full+1)*n),h(full+1),f(n);h[0]=1;
 for(unsigned s=1;s<=full;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q);unsigned previous=s^(1u<<v);int64_t c=previous?0:1;
  for(unsigned p=previous&in[v];p;p&=p-1)c=add(c,d[previous*n+countr_zero(p)]);
  d[s*n+v]=c;h[s]=add(h[s],c);}
 for(unsigned s=1;s<=full;s++)for(unsigned q=s;q;q&=q-1){int v=countr_zero(q);f[v]=add(f[v],mul(d[s*n+v],h[full^s]));}
 cout<<h[full];for(auto v:f)cout<<' '<<v;cout<<'\n';
 // Explicit checked arithmetic also covers our sparse order-21 targets.
}}
