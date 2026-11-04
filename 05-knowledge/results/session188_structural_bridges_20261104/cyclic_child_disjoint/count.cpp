#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
uint64_t paths(const vector<unsigned>&R){unsigned n=R.size();if(n>17)throw runtime_error("order bound");unsigned N=1u<<n;vector<uint64_t>d(size_t(N)*n);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;for(unsigned S=1;S<N;S++){unsigned bs=S;while(bs){unsigned b=bs&-bs;bs-=b;unsigned i=__builtin_ctz(b);auto z=d[size_t(S)*n+i];if(!z)continue;unsigned rest=R[i]&~S;while(rest){unsigned t=rest&-rest;rest-=t;d[size_t(S|t)*n+__builtin_ctz(t)]+=z;}}}uint64_t h=0;for(unsigned i=0;i<n;i++)h+=d[size_t(N-1)*n+i];return h;}
int main(){unsigned jobs;if(!(cin>>jobs))return 1;while(jobs--){unsigned n;cin>>n;vector<unsigned>R(n);for(auto &v:R)cin>>v;cout<<paths(R)<<'\n';}}
