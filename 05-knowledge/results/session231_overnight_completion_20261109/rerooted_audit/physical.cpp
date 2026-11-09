#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
using namespace std;
using U=uint64_t;
using Wide=__uint128_t;
U checked(Wide value){if(value>UINT64_MAX)throw runtime_error("uint64 overflow");return static_cast<U>(value);}
int main(){int n;while(cin>>n){vector<U>a(n);for(auto &x:a)cin>>x; if(n<1 || n>25)throw runtime_error("size guard"); U N=U(1)<<n,full=N-1;vector<uint32_t>d(N*n);vector<U>h(N);h[0]=1; vector<U>pred(n);for(int v=0;v<n;v++)for(int u=0;u<n;u++)if(a[u]>>v&1)pred[v]|=U(1)<<u;
for(U s=1;s<N;s++){U vs=s;while(vs){int v=__builtin_ctzll(vs);vs&=vs-1;U t=s^(U(1)<<v),z=0;if(!t)z=1;U ps=t&pred[v];while(ps){int u=__builtin_ctzll(ps);ps&=ps-1;z+=d[t*n+u];}if(z>UINT32_MAX)throw runtime_error("endpoint overflow");d[s*n+v]=z;h[s]=checked(Wide(h[s])+z);}}
vector<U>E(n),cov(n);for(int v=0;v<n;v++)E[v]=d[full*n+v];fill(d.begin(),d.end(),0);
for(U s=1;s<N;s++){U vs=s;while(vs){int v=__builtin_ctzll(vs);vs&=vs-1;U t=s^(U(1)<<v),z=0;if(!t)z=1;U ps=t&a[v];while(ps){int u=__builtin_ctzll(ps);ps&=ps-1;z+=d[t*n+u];}if(z>UINT32_MAX)throw runtime_error("endpoint overflow");d[s*n+v]=z;U product=checked(Wide(z)*h[full^s]);cov[v]=checked(Wide(cov[v])+product);}}
for(int v=0;v<n;v++)cout<<h[full]<<' '<<E[v]<<' '<<cov[v]<<' '<<h[full^(U(1)<<v)]<<'\n';cout<<"END\n";
}}
