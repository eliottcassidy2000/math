#include <iostream>
#include <vector>
#include <stdexcept>
using namespace std;
// Independent finite word automaton, not matrix inversion or fixed-multiplicity kernels.
int main(){int m,L,a,c;while(cin>>m>>L>>a>>c){if(m<1||m>8||L<1||L>12)throw runtime_error("bounds");int n=m+2,M=(1<<n)-1,H=m+L+1,mod=1<<L;vector<int> adj(n);for(int i=0;i<m;i++)cin>>adj[i];for(int i=0;i<m;i++){if(a>>i&1)adj[i]|=1<<m;else adj[m]|=1<<i;if(c>>i&1)adj[i]|=1<<(m+1);else adj[m+1]|=1<<i;}adj[m]|=1<<(m+1);
size_t width=(size_t)(M+1)*n*H;vector<int> dp((size_t)L*width);auto ix=[&](int e,int mask,int j,int h){return (size_t)e*width+((size_t)mask*n+j)*H+h;};for(int j=0;j<n;j++)dp[ix(0,1<<j,j,0)]=1;
for(int e=0;e<L;e++)for(int mask=1;mask<=M;mask++)for(int j=0;j<n;j++)if(mask>>j&1)for(int h=0;h<H;h++){int v=dp[ix(e,mask,j,h)];if(!v)continue;for(int k=0;k<n;k++){bool seen=mask>>k&1;if(seen&&k>=m)continue;int ne=e+seen;if(ne>=L)continue;int nm=mask|(1<<k);int w=seen?((-2*v)&(mod-1)):v;if((adj[j]>>k)&1){int&z=dp[ix(ne,nm,k,h)];z=(z+w)&(mod-1);}if(h+1<H){int&z=dp[ix(ne,nm,k,h+1)];z=(z+w)&(mod-1);}}}
for(int h=0;h<H;h++){int ans=0;for(int e=0;e<L;e++)for(int j=0;j<n;j++)ans=(ans+dp[ix(e,M,j,h)])&(mod-1);if(h)cout<<' ';cout<<ans;}cout<<'\n';}}
