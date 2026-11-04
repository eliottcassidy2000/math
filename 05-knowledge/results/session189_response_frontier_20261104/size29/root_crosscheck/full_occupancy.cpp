#include <iostream>
#include <vector>
#include <cstdint>
using namespace std;
// Independent mixed-radix full occupancy recurrence. No repeat bound or
// modularly weighted transitions: DP counts literal admissible run words.
// Child cover weights are applied only after all word states are counted.
int main(){int ncases;cin>>ncases;while(ncases--){int id,a,c;cin>>id>>a>>c;vector<int> q(9);for(int&x:q)cin>>x;vector<int> cap(11),stride(11),prof[11];for(int i=0;i<11;i++){cin>>cap[i];prof[i].resize(cap[i]+1);for(int j=1;j<=cap[i];j++)cin>>prof[i][j];}int total=1;for(int i=0;i<11;i++){stride[i]=total;total*=cap[i]+1;}cout<<id<<" "<<a<<" "<<c;for(int op=0;op<2;op++){vector<int> adj(11);for(int i=0;i<9;i++){adj[i]=op?(q[i]^511^(1<<i)):q[i];if(a>>i&1)adj[i]|=1<<9;if(c>>i&1)adj[i]|=1<<10;}adj[9]=(511^a)|(1<<10);adj[10]=511^c;vector<uint8_t> dp((size_t)total*11);for(int i=0;i<11;i++)dp[(size_t)stride[i]*11+i]=1;int result=0;for(int st=1;st<total;st++){int occ[11],available=0,weight=1;for(int i=0;i<11;i++){occ[i]=(st/stride[i])%(cap[i]+1);if(occ[i]<cap[i])available|=1<<i;if(!occ[i])weight=0;else weight=weight*prof[i][occ[i]]%8;}for(int end=0;end<11;end++){int val=dp[(size_t)st*11+end];if(!val)continue;result=(result+weight*val)%8;int mask=adj[end]&available;for(int next=0;next<11;next++)if(mask>>next&1){auto &cell=dp[(size_t)(st+stride[next])*11+next];cell=(cell+val)%8;}}}cout<<" "<<result;}cout<<"\n";}}
