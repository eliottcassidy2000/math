// Independent forward endpoint DP and two-bin joint port responses.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using I=int64_t;
int main(){int n,r,b;while(cin>>n>>r>>b){vector<unsigned>T(n);for(auto&x:T)cin>>x;int full=(1<<n)-1;
 vector<I>e((full+1)*n),s((full+1)*n),h(full+1);h[0]=1;
 for(int i=0;i<n;i++)e[(1<<i)*n+i]=s[(1<<i)*n+i]=1;
 for(int m=1;m<=full;m++)for(unsigned z=m;z;z&=z-1){int v=countr_zero(z);I a=e[m*n+v],c=s[m*n+v];
  h[m]+=a;for(unsigned q=T[v]&~m;q;q&=q-1)e[(m|(q&-q))*n+countr_zero(q)]+=a;
  for(unsigned q=(full^T[v]^(1<<v))&~m;q;q&=q-1)s[(m|(q&-q))*n+countr_zero(q)]+=c;}
 vector<int>p={r,b};I P=0;vector<I>E(2),S(2),W(4);
 for(int m=0;m<=full;m++){P+=h[m]*h[full^m];for(int i=0;i<2;i++){E[i]+=e[m*n+p[i]]*h[full^m];S[i]+=s[m*n+p[i]]*h[full^m];for(int j=0;j<2;j++)W[2*i+j]+=e[m*n+p[i]]*s[(full^m)*n+p[j]];}}
 cout<<h[full]<<' '<<P;for(int i=0;i<2;i++)cout<<' '<<e[full*n+p[i]]<<' '<<s[full*n+p[i]]<<' '<<h[full^(1<<p[i])]<<' '<<E[i]<<' '<<S[i];for(I x:W)cout<<' '<<x;cout<<'\n';
}}
