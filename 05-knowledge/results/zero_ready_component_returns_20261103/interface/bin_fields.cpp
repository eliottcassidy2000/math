// Exact small-component endpoint and three-bin fields, n<=17.
// Ordered-bin counts are bounded by binom(n+2,2)*n! < 2^63 at n<=17.
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std; using I=int64_t;
int main(){int n,w,u;while(cin>>n>>w>>u){
 if(n<2||n>17||w<0||w>=n||u<0||u>=n||u==w)throw runtime_error("invalid source or ports");
 int N=1<<n,full=N-1;vector<unsigned>T(n);for(auto&x:T)cin>>x;
 for(int i=0;i<n;i++){if(T[i]>>n||((T[i]>>i)&1))throw runtime_error("invalid row");
  for(int j=0;j<i;j++)if(((T[i]>>j)&1)+((T[j]>>i)&1)!=1)throw runtime_error("invalid pair");}
 vector<I> e(N*n),s(N*n),h(N),p2(N);h[0]=1;
 for(int v=0;v<n;v++)e[(1<<v)*n+v]=s[(1<<v)*n+v]=1;
 for(int mask=1;mask<N;mask++)for(unsigned vs=mask;vs;vs&=vs-1){int v=countr_zero(vs);
  I a=e[mask*n+v],b=s[mask*n+v];h[mask]+=a;
  for(unsigned z=T[v]&~mask;z;z&=z-1)e[(mask|(z&-z))*n+countr_zero(z)]+=a;
  for(unsigned z=(full^T[v]^(1<<v))&~mask;z;z&=z-1)s[(mask|(z&-z))*n+countr_zero(z)]+=b;}
 for(int mask=0;mask<N;mask++){int part=mask;while(true){p2[mask]+=h[part]*h[mask^part];if(!part)break;part=(part-1)&mask;}}
 vector<I> E2(n),S2(n),E3(n);
 for(int mask=1;mask<N;mask++)for(int v=0;v<n;v++){
  E2[v]+=e[mask*n+v]*h[full^mask];S2[v]+=s[mask*n+v]*h[full^mask];E3[v]+=e[mask*n+v]*p2[full^mask];}
 auto same_path=[&](int a,int b){vector<I>d(N*n);d[(1<<a)*n+a]=1;
  for(int mask=1;mask<N;mask++)for(unsigned vs=mask;vs;vs&=vs-1){int v=countr_zero(vs);
   for(unsigned z=T[v]&~mask;z;z&=z-1)d[(mask|(z&-z))*n+countr_zero(z)]+=d[mask*n+v];}
  I a2=0,a3=0;for(int mask=1;mask<N;mask++){a2+=d[mask*n+b]*h[full^mask];a3+=d[mask*n+b]*p2[full^mask];}
  return pair<I,I>(a2,a3);};
 auto awu=same_path(w,u),auw=same_path(u,w);
 I Bwu=0,G=0,D=0,SE=0,Fdelete=0;
 for(int mask=1;mask<N;mask++){
  int rest=full^mask;Bwu+=e[mask*n+u]*s[rest*n+w];D+=e[mask*n+w]*e[rest*n+u];
  if(!(mask>>u&1))Fdelete+=e[mask*n+w]*h[(full^(1<<u))^mask];
  if(!(mask>>u&1)||(mask>>w&1))continue;
  for(int part=rest;part;part=(part-1)&rest)if(part>>w&1){
   I factor=e[part*n+w]*h[rest^part];G+=e[mask*n+u]*factor;SE+=s[mask*n+u]*factor;}}
 vector<I>out={h[full],s[full*n+w],h[full^(1<<w)],E2[w],e[full*n+w],p2[full],S2[w],E3[w],p2[full^(1<<w)],
  e[full*n+u],s[full*n+u],h[full^(1<<u)],E2[u],S2[u],awu.first,auw.first,Bwu,G,D,auw.second,SE,Fdelete};
 for(size_t i=0;i<out.size();i++)cout<<(i?" ":"")<<out[i];cout<<'\n';
}}
