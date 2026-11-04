#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using U=uint64_t;
int main(){unsigned jobs;cin>>jobs;while(jobs--){unsigned n,h;cin>>n>>h;if(n>17||h<1||h>3)throw runtime_error("exact17-name first-three profile bound");vector<unsigned>Q(n);for(auto&x:Q)cin>>x;unsigned N=1u<<n;vector<U>d(size_t(N)*n*h);auto at=[&](unsigned S,unsigned v,unsigned k)->U&{return d[(size_t(S)*n+v)*h+k];};for(unsigned i=0;i<n;i++)at(1u<<i,i,0)=1;for(unsigned S=1;S<N;S++)for(unsigned last=0;last<n;last++)if(S>>last&1){unsigned js=(N-1)^S;while(js){unsigned b=js&-js;js-=b;unsigned j=__builtin_ctz(b);for(unsigned k=0;k<h;k++){U v=at(S,last,k);if(Q[last]>>j&1)at(S|b,j,k)+=v;if(k+1<h)at(S|b,j,k+1)+=v;}}}for(unsigned k=0;k<h;k++){U v=0;for(unsigned i=0;i<n;i++)v+=at(N-1,i,k);cout<<v<<' ';}cout<<'\n';}}
