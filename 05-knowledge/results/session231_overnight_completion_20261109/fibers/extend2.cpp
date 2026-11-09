#include <bits/stdc++.h>
using namespace std;using ll=long long;
struct Sig{ll a,g,s,e;int mask;bool operator<(const Sig& z)const{return tie(a,g,s,e)<tie(z.a,z.g,z.s,z.e);}bool eq(const Sig&z)const{return tie(a,g,s,e)==tie(z.a,z.g,z.s,z.e);}};
int main(){const int n=9,N=512,F=N-1;int base[]={64,121,35,37,13,81,28};long checked=0,unsat=0,dup=0,best=1;for(int code=0;code<(1<<15);code++){
 int adj[n]={};copy(base,base+7,adj);int z=0;for(int j=7;j<9;j++)for(int i=0;i<j;i++,z++){if(code>>z&1)adj[j]|=1<<i;else adj[i]|=1<<j;}
 ll en[N][n]={},st[N][n]={},h[N]={};h[0]=1;
 for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(s>>v&1){int r=s^(1<<v);if(!r)en[s][v]=st[s][v]=1;else for(int u=0;u<n;u++)if(r>>u&1){if(adj[u]>>v&1)en[s][v]+=en[r][u];if(adj[v]>>u&1)st[s][v]+=st[r][u];}h[s]+=en[s][v];}
 ll H=h[F],Q=0,aw[n]={},gw[n]={};for(int s=0;s<N;s++){Q+=h[s]*h[F^s];for(int v=0;v<n;v++){aw[v]+=st[s][v]*h[F^s];gw[v]+=en[s][v]*h[F^s];}}
 checked++;if(gcd(H,Q)!=1)continue;unsat++;Sig sig[N];sig[0]={0,0,0,0,0};for(int m=1;m<N;m++){int v=__builtin_ctz((unsigned)m),r=m^(1<<v);sig[m]={sig[r].a+aw[v],sig[r].g+gw[v],sig[r].s+st[F][v],sig[r].e+en[F][v],m};}sort(sig,sig+N);
 bool has=false;for(int i=1;i<N;i++)if(sig[i].eq(sig[i-1]))has=true;if(!has)continue;dup++;cerr<<"duplicate "<<code<<" H "<<H<<" Q "<<Q<<endl;
 ll K[n][n]={};for(int s=1;s<F;s++)for(int i=0;i<n;i++)if(s>>i&1)for(int j=0;j<n;j++)if(!(s>>j&1))K[i][j]+=en[s][i]*st[F^s][j];
 for(int i=0;i<N;){int j=i+1;while(j<N&&sig[j].eq(sig[i]))j++;if(j-i>1){map<ll,int> vals;for(int k=i;k<j;k++){ll t=0;int m=sig[k].mask;for(int u=0;u<n;u++)if(m>>u&1)for(int v=0;v<n;v++)if(m>>v&1)t+=K[u][v];vals[t]=m;}for(auto [t,m]:vals){int length=1;while(vals.count(t+2*length))length++;if(length>=2){best=max<long>(best,length);cout<<"{\"code\":"<<code<<",\"H\":"<<H<<",\"Q\":"<<Q<<",\"run\":"<<length<<",\"adj\":[";for(int k=0;k<n;k++)cout<<(k?",":"")<<adj[k];cout<<"],\"values\":[";for(int k=0;k<length;k++)cout<<(k?",":"")<<"["<<t+2*k<<","<<vals[t+2*k]<<"]";cout<<"]}"<<endl;}}}i=j;}
 }cerr<<"checked "<<checked<<" unsaturated "<<unsat<<" with_duplicate_interface "<<dup<<" best "<<best<<endl;
}
