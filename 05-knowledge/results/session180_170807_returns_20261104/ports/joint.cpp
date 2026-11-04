#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
using U=uint64_t;
int main(){int n,k;cin>>n>>k;vector<unsigned>T(n),ports(k);for(auto &v:T)cin>>v;for(auto &v:ports)cin>>v;unsigned all=(1u<<n)-1;size_t N=size_t(all+1)*n;vector<U>E(N),S(N),H(all+1);H[0]=1;for(int v=0;v<n;v++)E[(1u<<v)*n+v]=S[(1u<<v)*n+v]=1;
for(unsigned m=1;m<=all;m++){unsigned vs=m;while(vs){int v=__builtin_ctz(vs);vs&=vs-1;unsigned rem=m^(1u<<v);if(rem){unsigned q=rem&(~T[v]);U z=0;while(q){int u=__builtin_ctz(q);q&=q-1;z+=E[size_t(rem)*n+u];}E[size_t(m)*n+v]=z;q=rem&T[v];z=0;while(q){int u=__builtin_ctz(q);q&=q-1;z+=S[size_t(rem)*n+u];}S[size_t(m)*n+v]=z;}H[m]+=E[size_t(m)*n+v];}}
cout<<"H "<<H[all]<<'\n';for(unsigned p:ports){U f=0,g=0;for(unsigned m=1;m<=all;m++)if(m>>p&1){f+=E[size_t(m)*n+p]*H[all^m];g+=S[size_t(m)*n+p]*H[all^m];}cout<<"PORT "<<p<<' '<<E[size_t(all)*n+p]<<' '<<S[size_t(all)*n+p]<<' '<<H[all^(1u<<p)]<<' '<<f<<' '<<g<<'\n';}
for(unsigned p:ports)for(unsigned q:ports)if(p!=q){U j=0;for(unsigned m=1;m<all;m++)if((m>>p&1)&&!(m>>q&1))j+=E[size_t(m)*n+p]*S[size_t(all^m)*n+q];cout<<"J "<<p<<' '<<q<<' '<<j<<'\n';}}
