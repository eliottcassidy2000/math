#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <bit>
using U=uint32_t;
using A=std::vector<U>;

A source(const A& a,const A& b) {
  A t={ (U(1)<<(a.size()+1))-2 };
  for(U u:a)t.push_back(u<<1|((U(1)<<b.size())-1)<<(a.size()+1));
  for(U v:b)t.push_back(v<<(a.size()+1)|1);
  return t;
}
A attach(A s,const A& r,U w) {
  int n=s.size();
  U m=(U(1)<<r.size())-1;
  for(U &x:s)x|=m<<n;
  for(U y:r)s.push_back(y<<n);
  s[0]^=U(1)<<(n+w);
  s[n+w]^=1;
  return s;
}
bool strong(const A&a){
  U q=1;
  for(size_t j=0;j<a.size();j++)
    for(size_t i=0;i<a.size();i++)if(q>>i&1)q|=a[i];
  if(q!=(U(1)<<a.size())-1)return false;
  q=1;
  for(size_t j=0;j<a.size();j++)
    for(size_t i=0;i<a.size();i++)if(q&a[i])q|=U(1)<<i;
  return q==(U(1)<<a.size())-1;
}
U hp(const A&a) {
  int n=a.size();
  U f=(U(1)<<n)-1;
  A incoming(n,0);
  for (int i=0;i<n;i++)for(int j=0;j<n;j++)if((a[i]>>j)&1)incoming[j]|=U(1)<<i;
  // A valid directed path on a prescribed subset extends to a valid full
  // Hamilton path by deterministic vertex insertion and retains its
  // induced ordering. Thus each endpoint count is no larger than H.
  // Still check every narrowing rather than presuming the expected result.
  A p(size_t(f+1)*n,0);
  for (int i=0;i<n;i++)p[size_t(U(1)<<i)*n+i]=1;
  for(U s=1;s<=f;s++)if(s&(s-1))
    for(U x=s;x;x&=x-1){
      U v=std::countr_zero(x),t=s^(U(1)<<v);
      uint64_t z=0;
      for(U q=t&incoming[v];q;q&=q-1)z+=p[size_t(t)*n+std::countr_zero(q)];
      if(z>UINT32_MAX)throw std::overflow_error("narrowing");
      p[size_t(s)*n+v]=static_cast<U>(z);
    }
  uint64_t res=0;
  for(int i=0;i<n;i++)res+=p[size_t(f)*n+i];
  if(res>UINT32_MAX)throw std::overflow_error("total");
  return static_cast<U>(res);
}
int main() {
  const A a={154,124,169,192,44,137,53,82};
  const A b={100,57,18,21,1,28,62};
  A x=attach(source(a,b), {14,20,24,18,1}, 0);
  A y=attach(source(b,a), {26,60,57,48,32,1}, 1);
  auto h=hp(x),k=hp(y);
  if(h!=227125167||k!=227125165||!strong(x)||!strong(y))
    throw std::logic_error("incorrect outputs");
  std::cout<<x.size()<<": "<<h<<"; "<<y.size()<<": "<<k<<"; both strong\n";
}
