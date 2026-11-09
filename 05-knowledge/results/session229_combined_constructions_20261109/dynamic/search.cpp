#include <array>
#include <map>
#include <vector>
#include <iostream>
using namespace std;
struct W {int n,z,m;};
int main(){map<array<long long,6>,map<long long,W>> bank;int best=0;for(int n=1;n<=6;n++){int N=1<<n,edges=n*(n-1)/2;for(int z=0;z<(1<<edges);z++){int adj[6]={},k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(z>>k&1)adj[i]|=1<<j;else adj[j]|=1<<i;
long long en[64][6]={},st[64][6]={},h[64]={1};for(int s=1;s<N;s++){for(int v=0;v<n;v++)if(s>>v&1){int r=s^(1<<v);if(!r)en[s][v]=st[s][v]=1;else for(int u=0;u<n;u++)if(r>>u&1){if(adj[u]>>v&1)en[s][v]+=en[r][u];if(adj[v]>>u&1)st[s][v]+=st[r][u];}h[s]+=en[s][v];}}
long long q=0,aw[6]={},gw[6]={},tw[6][6]={};for(int s=0;s<N;s++){q+=h[s]*h[N-1-s];for(int i=0;i<n;i++){aw[i]+=st[s][i]*h[N-1-s];gw[i]+=en[s][i]*h[N-1-s];for(int j=0;j<n;j++)tw[i][j]+=en[s][i]*st[N-1-s][j];}}
for(int m=0;m<N;m++){long long A=0,G=0,T=0,S=0,E=0;for(int i=0;i<n;i++)if(m>>i&1){A+=aw[i];G+=gw[i];S+=st[N-1][i];E+=en[N-1][i];for(int j=0;j<n;j++)if(m>>j&1)T+=tw[i][j];}array<long long,6>key={h[N-1],q,A,G,S,E};auto &f=bank[key];f.emplace(T,W{n,z,m});int lo=T,hi=T;while(f.count(lo-2))lo-=2;while(f.count(hi+2))hi+=2;int run=(hi-lo)/2+1;if(run>best){best=run;cout<<"RECORD "<<n<<" "<<best<<" KEY";for(auto v:key)cout<<" "<<v;cout<<" WITNESSES";for(int t=lo;t<=hi;t+=2){auto w=f[t];cout<<" ("<<t<<","<<w.n<<","<<w.z<<","<<w.m<<")";}cout<<endl;}}}cerr<<"n="<<n<<" groups="<<bank.size()<<" best="<<best<<endl;}}
