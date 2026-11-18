#include <iostream>
#include <utility>
#include <vector>
using namespace std;
int h(const vector<int>& out){
 int n=out.size(), f=(1<<n)-1,p[256][8]{},r=0;
 for(int i=0;i<n;i++)p[1<<i][i]=1;
 for(int s=1;s<f;s++)for(int i=0;i<n;i++) if(p[s][i])
  for(int a=out[i]&~s;a;a&=a-1){int b=a&-a;p[s|b][__builtin_ctz(unsigned(b))]+=p[s][i];}
 for(int x:p[f])r+=x;
 return r;
}
vector<int> clone(const vector<int>& adj){
 int src[]={0,1,1,2,3,4,5,6};vector<int> res(8);
 for(int i=0;i<8;i++)for(int j=0;j<8;j++)if((src[i]==src[j]&&i<j)||(adj[src[i]]>>src[j]&1))res[i]|=1<<j;
 return res;
}
int main(){
 vector<pair<int,int>> es;
 for(int i=0;i<7;i++)for(int j=i+1;j<7;j++)es.push_back({i,j});
 int z=0,y=0;
 for(int t=0;t<(1<<es.size());t++){
  vector<int> adj(7);
  for(unsigned i=0;i<es.size();i++){auto [u,v]=es[i];if(t>>i&1)swap(u,v);adj[u]|=1<<v;}
  if(h(adj)!=35 || h(clone(adj))!=77) continue;
  ++z;
  vector<pair<int,int>> ed;
  for(int u=0;u<7;u++)if(adj[1]>>u&1)for(int v=0;v<7;v++)
   if((adj[v]>>1&1) && (adj[u]>>v&1)) ed.push_back({u,v});
  bool two=false;
  for(auto [a,b]:ed)for(auto [c,d]:ed)if(a!=c&&b!=d)two=true;
  if(!two && ++y==1){cout<<t<<": ";for(int a:adj)cout<<a<<' ';cout<<endl;}
 }
 cout<<z<<" interpolation hits, "<<y<<" with at most two visits"<<endl;
}
