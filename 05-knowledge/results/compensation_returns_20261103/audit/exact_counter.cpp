// Independent literal endpoint counter for compensation certificates.
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;using I=int64_t;
vector<I> paths(const vector<unsigned>&r){int n=r.size(),N=1<<n;vector<I>d(N*n);
 for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;
 for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(I z=d[m*n+v])
  for(unsigned q=r[v]&~m;q;q&=q-1){unsigned bit=q&-q;d[(m|bit)*n+countr_zero(bit)]+=z;}
 return d;
}
void emit(const vector<I>&v){cout<<'[';for(int i=0;i<(int)v.size();i++){if(i)cout<<',';cout<<v[i];}cout<<']';}
int main(int argc,char**argv){bool only_h=argc>1&&string(argv[1])=="h";int id,n;
 while(cin>>id>>n){if(n<1||n>19)throw runtime_error("declared exact order range1..19");
  vector<unsigned>r(n),rev(n);for(auto&v:r)cin>>v;unsigned full=(1u<<n)-1;
  for(int i=0;i<n;i++){if(r[i]&~full||r[i]>>i&1)throw runtime_error("row validity");
   for(int j=0;j<n;j++){if(i!=j&&((r[i]>>j&1)+(r[j]>>i&1))!=1)throw runtime_error("pair validity");if(r[j]>>i&1)rev[i]|=1u<<j;}}
  auto e=paths(r);int N=1<<n;vector<I>h(N);h[0]=1;
  for(int m=1;m<N;m++)for(int v=0;v<n;v++)h[m]+=e[m*n+v];
  cout<<"{\"id\":"<<id<<",\"H\":"<<h[full];
  if(!only_h){auto s=paths(rev);vector<I>E(n),S(n),F(n),G(n),Q(n);vector<vector<I>>J(n,vector<I>(n));I P=0;
   for(int m=0;m<N;m++){int q=full^m;P+=h[m]*h[q];
    for(int u=0;u<n;u++){F[u]+=e[m*n+u]*h[q];G[u]+=s[m*n+u]*h[q];
     if(e[m*n+u])for(int v=0;v<n;v++)J[u][v]+=e[m*n+u]*s[q*n+v];}}
   for(int v=0;v<n;v++){E[v]=e[full*n+v];S[v]=s[full*n+v];Q[v]=h[full^(1<<v)];}
   cout<<",\"P2\":"<<P<<",\"E1\":";emit(E);cout<<",\"S1\":";emit(S);cout<<",\"E2\":";emit(F);cout<<",\"S2\":";emit(G);cout<<",\"deletion_H\":";emit(Q);cout<<",\"J\":[";
   for(int u=0;u<n;u++){if(u)cout<<',';emit(J[u]);}cout<<']';
  }cout<<"}\n";
 }
}
