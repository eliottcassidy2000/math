#include <array>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using R=array<int,8>;
int main(int argc, char**) {
 map<int,int> want={{3,1},{5,1},{11,1},{13,5},{17,13},{19,7},{23,17},{29,11},{31,23}};
 if(argc>1)want={{5,3},{11,7},{13,9},{17,11},{19,13},{23,15},{29,19},{31,21}};
 for(int n=3;n<=7;n++){
  vector<pair<int,int>> free;
  for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)free.push_back({i,j});
  map<int,pair<int,int>> counts;
  for(int w=0;w<1<<free.size();w++){
   R r={};for(int i=0;i<n-1;i++)r[i]=1<<(i+1);
   for(int i=0;i<(int)free.size();i++){
     auto [a,b]=free[i];if(w>>i&1) r[a]|=1<<b; else r[b]|=1<<a;
   }
   int N=1<<n, a[256][8]={};
   for(int u=0;u<n;u++)a[1<<u][u]=1;
   for(int s=1;s<N-1;s++) for(int u=0;u<n;u++) if(a[s][u])
   for(int e=r[u]&~s;e;e&=e-1){int b=e&-e;a[s|b][__builtin_ctz(b)]+=a[s][u];}
   int h=0;for(int u=0;u<n;u++)h+=a[N-1][u];
   if(!want.count(h))continue;
   counts[h].first++;
   for(int s=1;s<N-1;s++){
    int k=0;for(int u=0;u<n;u++)k+=a[s][u];
    if(k==want[h]){counts[h].second++;break;}
   }
  }
  for(auto [p,v]:counts)cout<<n<<" "<<p<<" "<<want[p]<<" "<<v.first<<" "<<v.second<<"\n";
 }
}
