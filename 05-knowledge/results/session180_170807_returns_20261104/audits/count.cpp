#include <bits/stdc++.h>
using namespace std;
int main(){int n; while(cin>>n){ vector<unsigned> out(n),in(n);for(auto &x:out)cin>>x;for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(out[j]>>i&1)in[i]|=1u<<j; size_t N=1u<<n; vector<unsigned long long>d(N*n);for(int i=0;i<n;++i)d[(1u<<i)*n+i]=1;for(unsigned s=1;s<N;++s){unsigned v=s;while(v){int i=__builtin_ctz(v);v&=v-1;unsigned t=s^(1u<<i);if(!t)continue;unsigned p=t&in[i];while(p){int j=__builtin_ctz(p);p&=p-1;d[s*n+i]+=d[t*n+j];}}}unsigned long long h=0;for(int i=0;i<n;++i)h+=d[(N-1)*n+i];cout<<h<<'\n';}}
