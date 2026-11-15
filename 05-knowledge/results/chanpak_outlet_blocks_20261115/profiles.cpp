#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <utility>
using namespace std;
int main() {
 for(int n=2;n<=7;++n){
  set<pair<int,int>> profiles;
  int full=(1<<n)-1;
  for(int code=0;code<(1<<(n*(n-1)/2));++code){
   int r[7]{},irows[7]{};
   for(int i=0,k=0;i<n;++i)for(int j=i+1;j<n;++j,++k){
     if(code>>k&1)r[i]|=1<<j,irows[j]|=1<<i;
     else r[j]|=1<<i,irows[i]|=1<<j;
   }
   int seen=1;
   for(int old=0; old!=seen;){old=seen;for(int p=old;p;p&=p-1) seen|=r[__builtin_ctz(p)];}
   if(seen!=full) continue;
   seen=1;
   for(int old=0; old!=seen;){old=seen;for(int p=old;p;p&=p-1) seen|=irows[__builtin_ctz(p)];}
   if(seen!=full) continue;
   int ends[128][7]{}, h[128]{1};
   for(int mask=1;mask<=full;++mask){
    for(int p=mask;p;p&=p-1){
     int v=__builtin_ctz(p), old=mask^(1<<v);
     int s=old==0;
     for(int q=old&irows[v];q;q&=q-1)s+=ends[old][__builtin_ctz(q)];
     h[mask]+=(ends[mask][v]=s);
    }
   }
   int p=0;
   for(int m=0;m<=full;++m)p+=h[m]*h[full^m];
   profiles.insert({h[full],p});
  }
  cout<<n<<' '<<profiles.size()<<'\n';
  for(auto [h,p]:profiles) cout<<h<<' '<<p<<'\n';
 }
}
