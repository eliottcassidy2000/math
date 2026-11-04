#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
uint64_t count(const vector<unsigned>&T){unsigned n=T.size(),N=1u<<n;if(n>17)throw runtime_error("17! exact bound");vector<unsigned> pred(n);for(unsigned i=0;i<n;i++)for(unsigned j=0;j<n;j++)if(T[j]>>i&1)pred[i]|=1u<<j;vector<uint64_t>d(size_t(N)*n);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;for(unsigned s=1;s<N;s++){if(!(s&(s-1)))continue;unsigned ends=s;while(ends){unsigned bit=ends&-ends;ends-=bit;unsigned v=__builtin_ctz(bit),rest=s^bit,p=pred[v]&rest;while(p){unsigned b=p&-p;p-=b;d[size_t(s)*n+v]+=d[size_t(rest)*n+__builtin_ctz(b)];}}}uint64_t z=0;for(unsigned i=0;i<n;i++)z+=d[size_t(N-1)*n+i];return z;}
int main(){unsigned jobs;if(!(cin>>jobs))return 1;while(jobs--){unsigned n;cin>>n;vector<unsigned>T(n);for(auto&r:T)cin>>r;cout<<count(T)<<'\n';}}
