#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using U=uint64_t;
int main(){unsigned jobs;cin>>jobs;while(jobs--){unsigned n;cin>>n;if(n>17)throw runtime_error("17! binom(16,4) bound");vector<unsigned>T(n),pred(n);for(auto&r:T)cin>>r;for(unsigned i=0;i<n;i++)for(unsigned j=0;j<n;j++)if(T[j]>>i&1)pred[i]|=1u<<j;unsigned N=1u<<n;vector<U>d(size_t(N)*n),H(N);for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;for(unsigned S=1;S<N;S++){unsigned ends=S;while(ends){unsigned b=ends&-ends;ends-=b;unsigned i=__builtin_ctz(b),rest=S^b;if(rest){unsigned ps=pred[i]&rest;while(ps){unsigned p=ps&-ps;ps-=p;d[size_t(S)*n+i]+=d[size_t(rest)*n+__builtin_ctz(p)];}}H[S]+=d[size_t(S)*n+i];}}vector<U>prev(N);prev[0]=1;for(unsigned h=1;h<=5;h++){vector<U>next(N);for(unsigned S=1;S<N;S++)for(unsigned A=S;A;A=(A-1)&S)next[S]+=H[A]*prev[S^A];cout<<next[N-1]<<(h==5?'\n':' ');prev.swap(next);}}}
