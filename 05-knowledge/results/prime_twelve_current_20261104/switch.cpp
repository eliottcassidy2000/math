#include <iostream>
#include <vector>
#include <cstdint>
using I=std::uint64_t;
bool prime(const std::vector<unsigned>& a){
 unsigned n=a.size(),full=(1u<<n)-1;
 for(unsigned m=1;m<full;++m)if(__builtin_popcount(m)>1){
  bool yes=true;
  for(unsigned u=0;u<n;++u)if(!(m>>u&1)){
   unsigned out=a[u]&m;if(out && out!=m){yes=false;break;}
  }
  if(yes)return false;
 }
 return true;
}
bool strong(std::vector<unsigned> a){
 for(unsigned k=0;k<a.size();++k)for(unsigned u=0;u<a.size();++u)if(a[u]>>k&1)a[u]|=a[k];
 for(unsigned u=0;u<a.size();++u)if((a[u]|1u<<u)!=(1u<<a.size())-1)return false;
 return true;
}
int main(){
 const unsigned n=12,full=4095;std::vector<unsigned> base(n);
 const unsigned residues=(1<<1)|(1<<3)|(1<<4)|(1<<5)|(1<<9);
 for(unsigned u=0;u<11;++u)for(unsigned v=0;v<11;++v)if(residues>>((v+11-u)%11)&1)base[u]|=1<<v;
 base[11]=2047;
 for(unsigned flip=0;flip<2048;++flip){
  std::vector<unsigned>a=base;
  for(unsigned u=0;u<n;++u)a[u]^=(flip^((flip>>u&1)?full:0))&~(1u<<u);
  std::vector<I>s((full+1)*n),e(s.size());
  for(unsigned u=0;u<n;++u)s[(1u<<u)*n+u]=e[(1u<<u)*n+u]=1;
  for(unsigned mask=1;mask<=full;++mask)for(unsigned u=0;u<n;++u)if(mask>>u&1){
   unsigned rest=mask^(1<<u);if(!rest)continue;
   for(unsigned v=0;v<n;++v)if(rest>>v&1){
    if(a[u]>>v&1)s[mask*n+u]+=s[rest*n+v];
    else e[mask*n+u]+=e[rest*n+v];
   }
  }
  I h=0;unsigned sp=0,ep=0;
  for(unsigned u=0;u<n;++u){h+=s[full*n+u];sp|=(s[full*n+u]&1)<<u;ep|=(e[full*n+u]&1)<<u;}
  std::cout<<flip<<" "<<h<<" "<<sp<<" "<<ep<<" "<<prime(a)<<" "<<strong(a);
  for(unsigned u=0;u<n;++u)std::cout<<" "<<s[full*n+u]<<" "<<e[full*n+u];
  std::vector<I> profile(n+1);
  for(unsigned mask=1;mask<=full;++mask)for(unsigned u=0;u<n;++u)
   profile[__builtin_popcount(mask)]+=s[mask*n+u];
  for(unsigned k=1;k<=n;++k)std::cout<<" "<<profile[k];
  std::cout<<"\n";
 }
}

