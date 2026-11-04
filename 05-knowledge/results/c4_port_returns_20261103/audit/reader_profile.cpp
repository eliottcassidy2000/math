#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=int64_t;
vector<I> forward(const vector<unsigned>&t){int n=t.size(),N=1<<n;vector<I>d(N*n);
 for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(I c=d[m*n+v])
  for(unsigned q=t[v]&~m;q;q&=q-1){unsigned bit=q&-q;d[(m|bit)*n+countr_zero(bit)]+=c;}
 return d;
}
void emit_array(const vector<I>&a){cout<<'[';for(int i=0;i<(int)a.size();i++){if(i)cout<<',';cout<<a[i];}cout<<']';}
int main(){int id,n;while(cin>>id>>n){if(n<1||n>17)throw runtime_error("order bound");
 vector<unsigned>t(n),rev(n);for(auto&r:t)cin>>r;unsigned full=(1u<<n)-1;
 for(int i=0;i<n;i++){if(t[i]&~full||t[i]>>i&1)throw runtime_error("row validity");
  for(int j=0;j<n;j++){if(i!=j&&((t[i]>>j&1)+(t[j]>>i&1))!=1)throw runtime_error("pair validity");if(t[j]>>i&1)rev[i]|=1u<<j;}}
 auto e=forward(t),s=forward(rev);int N=1<<n;vector<I>h(N),E(n),S(n),F(n),G(n);h[0]=1;
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)h[m]+=e[m*n+v];
 vector<vector<I>>J(n,vector<I>(n));I P=0;
 for(int m=0;m<N;m++){int q=full^m;P+=h[m]*h[q];
  for(int u=0;u<n;u++){F[u]+=e[m*n+u]*h[q];G[u]+=s[m*n+u]*h[q];
   if(e[m*n+u])for(int v=0;v<n;v++)J[u][v]+=e[m*n+u]*s[q*n+v];}}
 for(int v=0;v<n;v++){E[v]=e[full*n+v];S[v]=s[full*n+v];}
 cout<<"{\"id\":"<<id<<",\"H\":"<<h[full]<<",\"P2\":"<<P<<",\"E1\":";emit_array(E);
 cout<<",\"S1\":";emit_array(S);cout<<",\"E2\":";emit_array(F);cout<<",\"S2\":";emit_array(G);cout<<",\"J\":[";
 for(int u=0;u<n;u++){if(u)cout<<',';emit_array(J[u]);}cout<<"]}\n";
}}
