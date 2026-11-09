#include <array>
#include <map>
#include <set>
#include <iostream>
#include <algorithm>
using namespace std;
int main(){map<array<int,6>,set<int>> fibers; for(int n=1;n<=6;n++){int m=1<<n,F=m-1; long long roots=0;for(int code=0;code<(1<<(n*(n-1)/2));code++) {int out[6]={},k=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(code>>k&1)out[i]|=1<<j;else out[j]|=1<<i;
int en[64][6]={},st[64][6]={},h[64]={};h[0]=1;
for(int U=1;U<m;U++)for(int v=0;v<n;v++)if(U>>v&1){int W=U^(1<<v);if(!W) en[U][v]=st[U][v]=1;else for(int w=0;w<n;w++)if(W>>w&1){if(out[w]>>v&1)en[U][v]+=en[W][w];if(out[v]>>w&1)st[U][v]+=st[W][w];}h[U]+=en[U][v];}
int q=0;for(int U=0;U<m;U++)q+=h[U]*h[F^U];
for(int C=0;C<m;C++){int A=0,G=0,T=0,S=0,E=0;for(int v=0;v<n;v++)if(C>>v&1){S+=st[F][v];E+=en[F][v];}for(int U=0;U<m;U++){int ep=0,sq=0;for(int v=0;v<n;v++)if(C>>v&1){ep+=en[U][v];sq+=st[F^U][v];}A+=h[U]*sq;G+=ep*h[F^U];T+=ep*sq;}fibers[{h[F],q,A,G,S,E}].insert(T);roots++;}}
int best=0;for(auto &kv:fibers)for(int x:kv.second){int run=1;while(kv.second.count(x+2*run))run++;best=max(best,run);}cout<<n<<" "<<roots<<" "<<fibers.size()<<" "<<best<<"\n";}}
