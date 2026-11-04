#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
using namespace std;
int main(){unsigned n;while(cin>>n){if(n>17)throw runtime_error("exact bound 17!");vector<unsigned>r(n);for(auto &v:r)cin>>v;unsigned N=1u<<n;vector<uint64_t>d(size_t(N)*n);for(unsigned v=0;v<n;v++)d[size_t(1u<<v)*n+v]=1;for(unsigned S=1;S<N;S++){for(unsigned v=0;v<n;v++){auto x=d[size_t(S)*n+v];if(!x)continue;unsigned next=r[v]&~S;while(next){unsigned b=next&-next;next-=b;d[size_t(S|b)*n+__builtin_ctz(b)]+=x;}}}uint64_t answer=0;for(unsigned v=0;v<n;v++)answer+=d[size_t(N-1)*n+v];cout<<answer<<'\n';}}
