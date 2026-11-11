// Literal endpoint DP for a named rooted tournament.
// mode0 uses the proved one-out-neighbour identity for U; mode1 counts U.
#include <iostream>
#include <vector>
using namespace std;using I=long long;
int main(){
 int n,r,mode;if(!(cin>>n>>r>>mode))return 1;vector<int>a(n);for(auto&x:a)cin>>x;
 int M=(1<<n)-1;vector<I>dp((size_t)(M+1)*n),h(M+1);
 for(int i=0;i<n;++i)dp[(size_t)(1<<i)*n+i]=1;
 for(int m=1;m<=M;++m)for(int i=0;i<n;++i)if(I x=dp[(size_t)m*n+i]){
  h[m]+=x;int v=a[i]&~m&M;
  while(v){int b=v&-v;v-=b;dp[(size_t)(m|b)*n+__builtin_ctz((unsigned)b)]+=x;}
 }
 h[0]=1;I H=h[M],E=dp[(size_t)M*n+r],d=h[M^(1<<r)],U=0;
 if(mode){
  fill(dp.begin(),dp.end(),0);dp[(size_t)(1<<r)*n+r]=1;
  for(int m=1;m<=M;++m)for(int i=0;i<n;++i)if(I x=dp[(size_t)m*n+i]){
   U+=x*h[M^m];int v=a[i]&~m&M;
   while(v){int b=v&-v;v-=b;dp[(size_t)(m|b)*n+__builtin_ctz((unsigned)b)]+=x;}
  }
 }else{if(__builtin_popcount((unsigned)a[r])!=1)return 2;U=H+d-E;}
 cout<<H<<' '<<E<<' '<<U<<' '<<d<<'\n';
}
