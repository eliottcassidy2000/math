#include <cstdint>
#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
uint64_t count(const vector<unsigned>& a){
 const unsigned n=a.size(),N=1u<<n,full=N-1;
 vector<unsigned> pred(n);for(unsigned v=0;v<n;v++)pred[v]=full^a[v]^(1u<<v);
 vector<uint64_t> d(size_t(N)*n);
 for(unsigned v=0;v<n;v++)d[size_t(1u<<v)*n+v]=1;
 for(unsigned mask=1;mask<N;mask++){
  unsigned ends=mask;
  while(ends){unsigned v=__builtin_ctz(ends);ends&=ends-1;unsigned rest=mask^(1u<<v);if(!rest)continue;unsigned p=rest&pred[v];uint64_t z=0;while(p){unsigned u=__builtin_ctz(p);p&=p-1;z+=d[size_t(rest)*n+u];}d[size_t(mask)*n+v]=z;}
 }
 uint64_t result=0;for(unsigned v=0;v<n;v++)result+=d[size_t(full)*n+v];return result;
}
int main(){unsigned cases;if(!(cin>>cases))return 1;while(cases--){unsigned n,a,c;cin>>n;if(n>15)throw runtime_error("unsupported old size");vector<unsigned>T(n+2);for(unsigned u=0;u<n;u++)cin>>T[u];cin>>a>>c;T[n]=1u<<(n+1);T[n+1]=0;for(unsigned u=0;u<n;u++){if(a>>u&1)T[u]|=1u<<n;else T[n]|=1u<<u;if(c>>u&1)T[u]|=1u<<(n+1);else T[n+1]|=1u<<u;}cout<<count(T)<<'\n';}}
