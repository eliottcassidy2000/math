// Independent full eight-field seed readout, with exact subset convolutions.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std; using I=int64_t;
int main(){int id,n,w;while(cin>>id>>n>>w){
 if(n<1||n>13||w<0||w>=n)return 2;
 vector<unsigned>r(n);for(auto&x:r)cin>>x;int N=1<<n,F=N-1;
 vector<I>e(N*n),s(N*n),h(N),p(N);h[0]=1;
 for(int M=1;M<N;M++)for(unsigned vs=M;vs;vs&=vs-1){
  int v=countr_zero(vs),R=M^(1<<v);I E=R?0:1,S=E;
  for(unsigned us=R;us;us&=us-1){int u=countr_zero(us);if(r[u]>>v&1)E+=e[R*n+u];if(r[v]>>u&1)S+=s[R*n+u];}
  e[M*n+v]=E;s[M*n+v]=S;h[M]+=E;
 }
 for(int M=0;M<N;M++){int U=M;while(true){p[M]+=h[U]*h[M^U];if(!U)break;U=(U-1)&M;}}
 I S2=0,E2=0,S3=0;
 for(int M=0;M<N;M++){S2+=s[M*n+w]*h[F^M];E2+=e[M*n+w]*h[F^M];S3+=s[M*n+w]*p[F^M];}
 cout<<"{\"id\":"<<id<<",\"fields\":["<<h[F]<<','<<e[F*n+w]<<','<<S2<<','<<E2<<','<<p[F^(1<<w)]<<','<<h[F^(1<<w)]<<','<<p[F]<<','<<S3<<"]}\n";
}}
