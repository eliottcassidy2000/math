// Independent literal-child forward path DP; does not use rooted F or h+F-b.
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
int main(){int n,b,id;long long count=0;
 while(cin>>n>>b>>id){vector<unsigned>t(n+1);for(int i=0;i<n;i++)cin>>t[i];
  t[0]|=1u<<n;t[n]=((1u<<n)-1)^1u;int N=n+1,full=(1<<N)-1;
  vector<int64_t>d((full+1)*N);
  for(int i=0;i<N;i++)d[(1<<i)*N+i]=1;
  for(int s=1;s<=full;s++)for(unsigned z=s;z;z&=z-1){int v=countr_zero(z);auto a=d[s*N+v];if(!a)continue;
   for(unsigned q=t[v]&~s;q;q&=q-1){int u=countr_zero(q);d[(s|(1<<u))*N+u]+=a;}}
  int64_t h=0;for(int i=0;i<N;i++)h+=d[full*N+i];
  if(h==6325)cout<<n<<' '<<b<<' '<<id<<'\n';count++;
 }
 cerr<<"literal children "<<count<<'\n';
}
