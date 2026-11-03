#include <bits/stdc++.h>
using namespace std;
int main(int argc,char**argv){if(argc!=2){cerr<<"usage: census n, 2<=n<=7\n";return 2;}int n=atoi(argv[1]);if(n<2||n>7){cerr<<"supported orders are 2..7\n";return 2;}int m=n*(n-1)/2;vector<pair<int,int>>es;int idx[10][10];for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){idx[i][j]=es.size();es.push_back({i,j});}
vector<vector<int>>maps;vector<int>flips;vector<int>p(n);iota(p.begin(),p.end(),0);do{vector<int>a(m);int f=0;for(int k=0;k<m;k++){auto[i,j]=es[k];int x=p[i],y=p[j];a[k]=idx[min(x,y)][max(x,y)];if(x>y)f|=1<<a[k];}maps.push_back(a);flips.push_back(f);}while(next_permutation(p.begin(),p.end()));
vector<int>cls(1<<m,-1),rep;for(int t=0;t<(1<<m);t++){if(cls[t]>=0)continue;int c=rep.size();rep.push_back(t);for(int k=0;k<maps.size();k++){int u=flips[k],s=t;while(s){int b=__builtin_ctz(s);u^=1<<maps[k][b];s&=s-1;}cls[u]=c;}}
vector<vector<int>>adj(rep.size());for(int c=0;c<rep.size();c++){for(int k=0;k<m;k++){int d=cls[rep[c]^(1<<k)];if(c!=d)adj[c].push_back(d);}sort(adj[c].begin(),adj[c].end());adj[c].erase(unique(adj[c].begin(),adj[c].end()),adj[c].end());}
vector<int>st;for(int c=0;c<rep.size();c++){vector<int>reach(n);for(int i=0;i<n;i++)reach[i]=1<<i;for(int k=0;k<m;k++){auto[i,j]=es[k];if(rep[c]>>k&1)reach[i]|=1<<j;else reach[j]|=1<<i;}for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];if(reach==vector<int>(n,(1<<n)-1))st.push_back(c);}
cout<<"{\"n\":"<<n<<",\"representatives\":[";for(int i=0;i<rep.size();i++){if(i)cout<<",";cout<<rep[i];}cout<<"],\"strong\":[";for(int i=0;i<st.size();i++){if(i)cout<<",";cout<<st[i];}cout<<"],\"edges\":[";bool first=true;for(int i=0;i<rep.size();i++)for(auto j:adj[i])if(i<j){if(!first)cout<<",";first=false;cout<<"["<<i<<","<<j<<"]";}cout<<"]}\n";
}
