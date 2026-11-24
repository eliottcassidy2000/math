#include <cstdint>
#include <cstdio>
#include <vector>
#include <utility>
using namespace std;
uint64_t h(const vector<int>& r, int cycle) {
  const int n=r.size(), full=(1<<n)-1;
  vector<uint64_t> F((full+1)*n);
  for(int i=0;i<n;i++) if(!cycle||!i) F[(1<<i)*n+i]=1;
  for(int S=1;S<full;S++){
    if(cycle && !(S&1))continue;
    for(int L=S;L;L&=L-1){
      int u=__builtin_ctz(unsigned(L));
      auto v=F[S*n+u];
      if(!v) continue;
      for(int a=r[u]&~S;a;a&=a-1)F[(S|(a&-a))*n+__builtin_ctz(unsigned(a))]+=v;
    }
  }
  uint64_t ans=0;
  for(int i=0;i<n;i++) if(!cycle || (r[i]&1)) ans+=F[full*n+i];
  return ans;
}
int main(){
  vector<pair<int,int>> V={{0,2},{0,3},{1,2},{1,3},{0,4},{1,4},{2,4},{3,4},{0,1},{2,3}};
  vector<int> rows(20);
  int n=V.size();
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
    rows[i]|=1<<j;rows[n+j]|=1<<(n+i);
  }
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){
    auto [a,b]=V[i];auto [c,d]=V[j];
    if(i==j||(a!=c&&a!=d&&b!=c&&b!=d))rows[i]|=1<<(n+j);else rows[n+j]|=1<<i;
  }
  for(int pass=0;pass<2;pass++){
    printf("Hamilton paths %llu, cycles %llu\n",(unsigned long long)h(rows,0),(unsigned long long)h(rows,1));
    vector<int> spine;for(int i=0;i<n;i++)spine.push_back(i);
    for(int i=2*n-1;i>=n;i--)spine.push_back(i);
    for(int i=0;i+1<2*n;i++){
      int u=spine[i],v=spine[i+1];rows[u]^=1<<v;rows[v]^=1<<u;
    }
  }
}
