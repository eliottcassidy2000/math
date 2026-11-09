#include <bits/stdc++.h>
using namespace std;using ll=long long;using U=unsigned long long;
int main(int argc,char**argv){int n=stoi(argv[1]),N=1<<n,F=N-1;map<array<ll,5>,map<ll,pair<U,int>>> fibers;map<array<ll,7>,pair<U,int>> selected;ifstream input(argv[2]);string line;long classes=0;
while(getline(input,line)){if(line.empty()||line[0]=='#')continue;istringstream is(line);U code;is>>code;if(!is)continue;classes++;int adj[n]={},z=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,z++){if(code>>z&1)adj[i]|=1<<j;else adj[j]|=1<<i;}
ll en[N][n]={},st[N][n]={},h[N]={};h[0]=1;for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(s>>v&1){int r=s^(1<<v);if(!r)en[s][v]=st[s][v]=1;else for(int u=0;u<n;u++)if(r>>u&1){if(adj[u]>>v&1)en[s][v]+=en[r][u];if(adj[v]>>u&1)st[s][v]+=st[r][u];}h[s]+=en[s][v];}
ll H=h[F],Q=0,aw[n]={},gw[n]={},K[n][n]={};for(int s=0;s<N;s++){Q+=h[s]*h[F^s];for(int i=0;i<n;i++){aw[i]+=st[s][i]*h[F^s];gw[i]+=en[s][i]*h[F^s];for(int j=0;j<n;j++)K[i][j]+=en[s][i]*st[F^s][j];}}
if(H==1||Q%H)continue;for(int m=0;m<N;m++){ll A=0,G=0,T=0,S=0,E=0;for(int i=0;i<n;i++)if(m>>i&1){A+=aw[i];G+=gw[i];S+=st[F][i];E+=en[F][i];for(int j=0;j<n;j++)if(m>>j&1)T+=K[i][j];}selected.emplace(array<ll,7>{H,Q,A,G,T,S,E},make_pair(code,m));}
}
for(auto &[s,w]:selected){cout<<"[";for(int i=0;i<7;i++)cout<<(i?",":"")<<s[i];cout<<","<<w.first<<","<<w.second<<"]\n";}cerr<<"classes "<<classes<<" states "<<selected.size()<<endl;
}
