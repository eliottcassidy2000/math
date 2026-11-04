#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
using I=std::int64_t;
std::vector<I> current(const std::vector<unsigned>& a){
 unsigned n=a.size(),end=(1u<<n)-1;
 std::vector<I> e((end+1)*n),s(e.size()),h(end+1);
 h[0]=1;
 for(unsigned u=0;u<n;++u)e[(1u<<u)*n+u]=s[(1u<<u)*n+u]=1;
 for(unsigned mask=1;mask<=end;++mask){
  for(unsigned u=0;u<n;++u)if(mask>>u&1){
   unsigned rest=mask^(1u<<u);
   if(rest)for(unsigned v=0;v<n;++v)if(rest>>v&1){
    if(a[v]>>u&1)e[mask*n+u]+=e[rest*n+v];
    else s[mask*n+u]+=s[rest*n+v];
   }
   h[mask]+=s[mask*n+u];
  }
 }
 std::vector<I> g(n);
 for(unsigned u=0;u<n;++u){
  g[u]=2*(e[end*n+u]-s[end*n+u]);
  for(unsigned mask=1;mask<end;++mask)if(mask>>u&1)
   g[u]+=(e[mask*n+u]-s[mask*n+u])*h[end^mask];
 }
 return g;
}
bool strong(std::vector<unsigned> a){
 unsigned n=a.size();
 for(unsigned k=0;k<n;++k)for(unsigned u=0;u<n;++u)if(a[u]>>k&1)a[u]|=a[k];
 for(unsigned u=0;u<n;++u)if((a[u]|(1u<<u))!=(1u<<n)-1)return false;
 return true;
}
bool prime(const std::vector<unsigned>& a){
 unsigned n=a.size(),end=(1u<<n)-1;
 for(unsigned m=1;m<end;++m)if(__builtin_popcount(m)>1){
  bool module=true;
  for(unsigned u=0;u<n;++u)if(!(m>>u&1)){
   auto x=a[u]&m;if(x && x!=m){module=false;break;}
  }
  if(module)return false;
 }
 return true;
}
int main(){
 unsigned n;std::uint64_t code;
 while(std::cin>>n>>code){
  std::vector<unsigned>a(n);
  unsigned bit=0;
  for(unsigned u=0;u<n;++u)for(unsigned v=u+1;v<n;++v,++bit)
   if(code>>bit&1)a[u]|=1u<<v;else a[v]|=1u<<u;
  auto g=current(a);
  if(std::any_of(g.begin(),g.end(),[](I x){return x!=0;}))continue;
  bool reg=true;
  for(unsigned u=0;u<n;++u)if(2*__builtin_popcount(a[u])!=n-1)reg=false;
  bool cards=true;
  for(unsigned gone=0;gone<n && cards;++gone){
   std::vector<unsigned>b(n-1);unsigned i=0;
   for(unsigned u=0;u<n;++u)if(u!=gone){
    unsigned j=0;for(unsigned v=0;v<n;++v)if(v!=gone){if(a[u]>>v&1)b[i]|=1u<<j;++j;}++i;
   }
   auto q=current(b);for(auto x:q)if(x)cards=false;
  }
  std::cout<<n<<" "<<code<<" "<<strong(a)<<" "<<reg<<" "<<prime(a)<<" "<<cards<<"\n";
 }
}

