#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
using namespace std;
int main(){int n;while(cin>>n){if(n<1||n>16)throw runtime_error("order bound");vector<int>out(n),in(n);for(int&i:out)cin>>i;for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(out[i]>>j&1)in[j]|=1<<i;int full=(1<<n)-1;vector<uint64_t>d((size_t)(full+1)*n);for(int mask=1;mask<=full;mask++){int bits=mask;while(bits){int j=__builtin_ctz(bits);bits&=bits-1;int prev=mask^(1<<j);uint64_t v=0;if(!prev)v=1;else{int x=in[j]&prev;while(x){int i=__builtin_ctz(x);x&=x-1;v+=d[(size_t)prev*n+i];}}d[(size_t)mask*n+j]=v;}}uint64_t ans=0;for(int i=0;i<n;i++)ans+=d[(size_t)full*n+i];cout<<ans<<'\n';}}
