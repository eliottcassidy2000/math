#include <bits/stdc++.h>
using namespace std;
int main(){
 int n;cin>>n;int F=(1<<n)-1;
 vector<int>A(n),in(n);for(int&i:A)cin>>i;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(A[i]>>j&1)in[j]|=1<<i;
 vector<uint32_t> e(size_t(n)*(1<<n));vector<uint32_t> h(1<<n);h[0]=1;
 for(int S=1;S<=F;S++)for(int r=S;r;r&=r-1){
  int j=__builtin_ctz(r),U=S^(1<<j);
  uint64_t u=0;
  if(!U)u=1;else for(int t=U&in[j];t;t&=t-1)u+=e[size_t(U)*n+__builtin_ctz(t)];
  if(u>UINT32_MAX || u+h[S]>UINT32_MAX) throw overflow_error("endpoint count"); e[size_t(S)*n+j]=u;h[S]+=u;
 }
 unsigned __int128 q=0,root=0;for(int S=0;S<=F;S++){q+=uint64_t(h[S])*h[F^S];if(S&1)root+=uint64_t(e[size_t(S)*n])*h[F^S];}
 if(q>UINT64_MAX || root>UINT64_MAX) throw overflow_error("cover count"); cout<<h[F]<<" "<<uint64_t(q)<<" "<<uint64_t(root)<<endl;
}
