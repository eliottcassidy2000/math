// Exhaustive retained-root one-birth responses for one named source.
// Input: n root followed by n outgoing-neighbour masks.
// Output: birth-cut H E U d for every cut.
#include <iostream>
#include <vector>
using namespace std;
using I=long long;
int main(){
 int n,r;if(!(cin>>n>>r))return 1;vector<int> old(n);for(auto &x:old)cin>>x;
 int N=n+1,M=(1<<N)-1;
 vector<I> dp((M+1)*N),sp((M+1)*N),h(M+1),starts(M+1);
 for(int cut=0;cut<(1<<n);++cut){
  vector<int> rows(N);for(int i=0;i<n;++i)rows[i]=old[i]|((cut&(1<<i))?0:(1<<n));rows[n]=cut;
  fill(dp.begin(),dp.end(),0);fill(sp.begin(),sp.end(),0);fill(h.begin(),h.end(),0);fill(starts.begin(),starts.end(),0);
  for(int i=0;i<N;++i)dp[(1<<i)*N+i]=1;sp[(1<<r)*N+r]=1;
  for(int m=1;m<=M;++m){
   for(int i=0;i<N;++i){
    I x=dp[m*N+i],y=sp[m*N+i];h[m]+=x;starts[m]+=y;
    if(!x&&!y)continue;
    int a=rows[i]&~m&M;
    while(a){int b=a&-a;a-=b;int j=__builtin_ctz((unsigned)b);dp[(m|b)*N+j]+=x;sp[(m|b)*N+j]+=y;}
   }
  }
  h[0]=1;I U=0;for(int m=1;m<=M;++m)U+=starts[m]*h[M^m];
  cout<<cut<<' '<<h[M]<<' '<<dp[M*N+r]<<' '<<U<<' '<<h[M^(1<<r)]<<'\n';
 }
}
