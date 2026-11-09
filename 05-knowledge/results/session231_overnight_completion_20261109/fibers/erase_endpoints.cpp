#include <bits/stdc++.h>
using namespace std;using ll=long long;
int main(){const int n=6,N=64,F=63;map<array<ll,5>,map<ll,pair<int,int>>> fibers;map<array<ll,7>,pair<int,int>> selected;
for(int code=0;code<32768;code++){int adj[n]={},z=0;for(int j=1;j<n;j++)for(int i=0;i<j;i++,z++){if(code>>z&1)adj[j]|=1<<i;else adj[i]|=1<<j;}
ll en[N][n]={},st[N][n]={},h[N]={};h[0]=1;for(int s=1;s<N;s++)for(int v=0;v<n;v++)if(s>>v&1){int r=s^(1<<v);if(!r)en[s][v]=st[s][v]=1;else for(int u=0;u<n;u++)if(r>>u&1){if(adj[u]>>v&1)en[s][v]+=en[r][u];if(adj[v]>>u&1)st[s][v]+=st[r][u];}h[s]+=en[s][v];}
ll H=h[F],Q=0,aw[n]={},gw[n]={},K[n][n]={};for(int s=0;s<N;s++){Q+=h[s]*h[F^s];for(int i=0;i<n;i++){aw[i]+=st[s][i]*h[F^s];gw[i]+=en[s][i]*h[F^s];for(int j=0;j<n;j++)K[i][j]+=en[s][i]*st[F^s][j];}}
for(int m=0;m<N;m++){ll A=0,G=0,T=0,S=0,E=0;for(int i=0;i<n;i++)if(m>>i&1){A+=aw[i];G+=gw[i];S+=st[F][i];E+=en[F][i];for(int j=0;j<n;j++)if(m>>j&1)T+=K[i][j];}if((H==29&&Q==320)||(H==45&&Q==432))selected.emplace(array<ll,7>{H,Q,A,G,T,S,E},make_pair(code,m));for(int v=0;v<2;v++)fibers[{H,Q,A,G,v}].emplace(T+S+v*E,make_pair(code,m));}
}
ofstream out("session231-dense-fibers/CORE_BANK.jsonl");for(auto &[s,w]:selected){out<<"[";for(int i=0;i<7;i++)out<<(i?",":"")<<s[i];out<<","<<w.first<<","<<w.second<<"]\n";}int best=0,unsatbest=0;for(auto &[key,vals]:fibers)for(auto [t,w]:vals){if(vals.count(t-2))continue;int length=1;while(vals.count(t+2*length))length++;bool unsat=gcd(key[0],key[1])==1;if(length>best||(unsat&&length>unsatbest)){best=max(best,length);if(unsat)unsatbest=max(unsatbest,length);cout<<"{\"key\":[";for(int i=0;i<5;i++)cout<<(i?",":"")<<key[i];cout<<"],\"run\":"<<length<<",\"unsaturated\":"<<unsat<<",\"witnesses\":[";for(int k=0;k<length;k++){auto [c,m]=vals[t+2*k];cout<<(k?",":"")<<"["<<t+2*k<<","<<c<<","<<m<<"]";}cout<<"]}"<<endl;}}
array<ll,5> key={45,432,129,129,0};cerr<<"target:";if(auto it=fibers.find(key);it!=fibers.end())for(auto [t,w]:it->second)cerr<<" "<<t;cerr<<"\nfibers "<<fibers.size()<<" best "<<best<<" unsatbest "<<unsatbest<<endl;
}
