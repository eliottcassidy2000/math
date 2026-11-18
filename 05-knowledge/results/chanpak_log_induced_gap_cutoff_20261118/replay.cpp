#include <iostream>
#include <map>
#include <utility>
#include <vector>
using namespace std;
int n, adj[8];
map<pair<int,int>,vector<int>> witnesses;
void record(int a,int b,int s) {
    if (b-a<2 || b-a>6 || witnesses.count({a,b})) return;
    vector<int> v={n,s};
    v.insert(v.end(),adj,adj+n);
    witnesses[{a,b}]=v;
}
bool strong(){
    for(int v=0;v<n;v++){
        int reach=1<<v,next;
        do {next=reach;for(int j=0;j<n;j++) if(reach>>j&1)next|=adj[j];
            if(next==reach) break;
            reach=next;}while(true);
        if(reach!=(1<<n)-1)return false;
    }return true;
}
int cycles;
void cycle_dfs(int v,int used,int len){
    if(len>=3 && len%2 && (adj[v]>>(n-1)&1))++cycles;
    for(int avail=adj[v]&~used;avail && cycles<=3;avail&=avail-1) {
        int w=__builtin_ctz(unsigned(avail));
        cycle_dfs(w,used|(1<<w),len+1);
    }
}
void audit(){
    int full=(1<<n)-1;
    if(n>=7){
        cycles=0;
        cycle_dfs(n-1,1<<(n-1),1);
        if(cycles!=3)return;
    }
    if(!strong())return;
    int dp[256][8]={}, h[256]={1};
    for(int s=1;s<=full;s++){
        for(int v=0;v<n;v++)if(s>>v&1){
            int t=s^(1<<v);
            if(!t)dp[s][v]=1;
            else for(int u=0;u<n;u++) if(adj[u]>>v&1)dp[s][v]+=dp[t][u];
            h[s]+=dp[s][v];
        }
    }
    if(n>=7)record(h[full^(1<<(n-1))],h[full],full^(1<<(n-1)));
    else for(int s=0;s<full;s++)record(h[s],h[full],s);
}
int main(){
    for(n=3;n<=8;n++){
        vector<pair<int,int>> edges;
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
            if(!(n>=7 && j==i+1 && j<n-1))edges.push_back({i,j});
        for(int b=0;b<(1<<edges.size());b++){
            for(int i=0;i<n;i++)adj[i]=(n>=7 && i+1<n-1)?1<<(i+1):0;
            for(unsigned k=0;k<edges.size();k++){
                auto [i,j]=edges[k];
                if(b>>k&1)adj[j]|=1<<i;else adj[i]|=1<<j;
            }
            audit();
        }
        cout<<"finished strong hosts of order "<<n<<'\n';
    }
    for(auto [p,v]:witnesses){
        cout<<p.first<<" -> "<<p.second<<": ";
        for(int x:v)cout<<x<<' ';
        cout<<'\n';
    }
}
