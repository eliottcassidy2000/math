#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
uint64_t count(const vector<unsigned>&G){unsigned n=G.size();if(n>16)throw runtime_error("exact uint64 order boundary");unsigned N=1u<<n;vector<unsigned>pred(n);for(unsigned i=0;i<n;i++)for(unsigned j=0;j<n;j++)if(G[i]>>j&1)pred[j]|=1u<<i;vector<uint64_t>D(size_t(N)*n);
for(unsigned M=1;M<N;M++){unsigned bs=M;while(bs){unsigned b=bs&-bs;bs-=b;unsigned v=__builtin_ctz(b),R=M^b;uint64_t value=0;if(!R)value=1;else{unsigned p=pred[v]&R;while(p){unsigned t=p&-p;p-=t;value+=D[size_t(R)*n+__builtin_ctz(t)];}}D[size_t(M)*n+v]=value;}}
uint64_t H=0;for(unsigned i=0;i<n;i++)H+=D[size_t(N-1)*n+i];return H;}
int main(){unsigned jobs;if(!(cin>>jobs))return 1;while(jobs--){unsigned n;cin>>n;vector<unsigned>G(n);for(auto&v:G)cin>>v;cout<<count(G)<<'\n';}}
