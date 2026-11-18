#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
using namespace std;
int main(){
 constexpr int n=6,L=1<<15,f=(1<<n)-1;
 vector<pair<int,int>> edges;
 for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)edges.push_back({u,v});
 vector<int> H(L);
 for(int z=0;z<L;z++){
  int m[n]{},d[64][n]{};
  for(int k=0;k<15;k++){
   auto [u,v]=edges[k];if(z>>k&1)swap(u,v);m[u]|=1<<v;
  }
  for(int u=0;u<n;u++)d[1<<u][u]=1;
  for(int s=1;s<f;s++)for(int u=0;u<n;u++)if(d[s][u])
   for(int t=m[u]&~s;t;t&=t-1)d[s|(t&-t)][__builtin_ctz(unsigned(t))]+=d[s][u];
  for(int u=0;u<n;u++)H[z]+=d[f][u];
 }

 for(int a=9;a<=19;a+=2)for(int b=a+2;b<=19;b+=2){
  int best=15,lo=-1,hi=-1;
  vector<int> C;
  for(int z=0;z<L;z++)if(H[z]==b)C.push_back(z);
  for(int z=0;z<L;z++)if(H[z]==a)
   for(int w:C){int d=__builtin_popcount(unsigned(w^z));if(d<best){best=d;lo=z;hi=w;if(best==1)break;}}
  cout<<a<<" "<<b<<" "<<best<<" "<<lo<<" "<<hi<<'\n';
 }
}
