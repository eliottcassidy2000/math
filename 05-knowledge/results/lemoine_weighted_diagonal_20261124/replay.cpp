#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
int main(int argc, char**argv) {
  const int top=argc>1?std::stoi(argv[1]):100000000;
  assert(top>=21 && top<=500000000);
  std::vector<bool> prime(top+1,true);
  prime[0]=prime[1]=false;
  for(int p=2; int64_t(p)*p<=top; ++p) if(prime[p])
    for(int64_t j=int64_t(p)*p;j<=top;j+=p) prime[j]=false;
  std::vector<int> ps;
  for(int p=2;p<=top;++p) if(prime[p]) ps.push_back(p);
  int uq=0, up=0, bal=0, no_bridge=0, no_offdiag=0, missing=0;
  for(int n=7;n<=top;n+=2) {
    int q_min=0, p_min=0, q_off=0, both=0;
    for(int q:ps) {
      if(2*q>=n) break;
      const int p=n-2*q;
      if(!prime[p]) continue;
      if(!q_min) q_min=q;
      if(q!=2 && q!=p){ q_off=q; break; }
    }
    for(int p:ps) {
      if(p>=n) break;
      if(p%2 && (n-p)%2==0 && prime[(n-p)/2]){ p_min=p; break; }
    }
    // Search only up to the first attained min(p,q) from either side.
    // For m increasing, test m as one endpoint and the other endpoint >= m.
    for(int m:ps) {
      if(3*m>n) break;
      if((m%2 && (n-m)%2==0 && prime[(n-m)/2]) ||
         prime[n-2*m]) { both=m; break; }
    }
    assert(both==std::min(p_min,q_min));
    if(!q_min) { ++missing; std::cout<<"missing "<<n<<"\n"; }
    if(!q_off) { ++no_offdiag; std::cout<<"no odd off-diagonal "<<n<<"\n"; }
    if(!q_off && (n%3!=0 || !prime[n/3])){
      ++no_bridge; std::cout<<"no odd-prime bridge "<<n<<"\n";
    }
    if(q_min>uq){ uq=q_min; std::cout<<"q-record "<<n<<" "<<uq<<" witness "<<n-2*uq<<"\n"; }
    if(p_min>up){ up=p_min; std::cout<<"p-record "<<n<<" "<<up<<" witness "<<(n-up)/2<<"\n"; }
    if(both>bal){ bal=both; std::cout<<"both-record "<<n<<" "<<bal<<"\n"; }
    if(n==7 || n==21)
      for(int q:ps) {
        if(2*q>=n) break;
        if(prime[n-2*q]) std::cout<<"pair "<<n<<" "<<n-2*q<<" "<<q<<"\n";
      }
  }
  std::cout<<"through "<<top<<"; missing="<<missing<<"; no odd bridge="<<no_bridge
           <<"; no odd off-diagonal="<<no_offdiag<<"\n";
  assert(!missing);
}
