#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <numeric>
using namespace std;
vector<int> outs(int k,int mask){vector<int>o(k);int t=0;for(int u=0;u<k;u++)for(int v=u+1;v<k;v++,t++){if(mask>>t&1)o[u]|=1<<v;else o[v]|=1<<u;}return o;}
long long resp(vector<int>o,int a,int b){int k=o.size(),n=k+2; o.resize(n);o[k]|=1<<(k+1);for(int u=0;u<k;u++){if(a>>u&1)o[u]|=1<<k;else o[k]|=1<<u;if(b>>u&1)o[u]|=1<<(k+1);else o[k+1]|=1<<u;}vector<long long>d((1<<n)*n);for(int u=0;u<n;u++)d[(1<<u)*n+u]=1;for(int s=1;s<(1<<n);s++)for(int u=0;u<n;u++)if(d[s*n+u]){int z=o[u]&~s;while(z){int v=__builtin_ctz(z);z&=z-1;d[(s|1<<v)*n+v]+=d[s*n+u];}}long long h=0;for(int u=0;u<n;u++)h+=d[((1<<n)-1)*n+u];return h;}
int main(int argc,char**argv){int k=argc>1?atoi(argv[1]):5,ne=k*(k-1)/2;vector<vector<int>>idx(k,vector<int>(k));int t=0;for(int u=0;u<k;u++)for(int v=u+1;v<k;v++)idx[u][v]=idx[v][u]=t++;map<vector<int>,vector<int>>gr;
for(int q=0;q<(1<<ne);q++){auto o=outs(k,q);vector<int>p(k);iota(p.begin(),p.end(),0);vector<int>s(1+k+ne),st(k),en(k);vector<vector<int>>c(k,vector<int>(k));do{int bad=0,badpos=-1;for(int j=1;j<k;j++)if(!(o[p[j-1]]>>p[j]&1)){bad++;badpos=j;}if(bad==0){s[0]++;st[p[0]]++;en[p[k-1]]++;}if(bad<=1)for(int j=1;j<k;j++)if(bad==0||j==badpos)c[p[j-1]][p[j]]++;}while(next_permutation(p.begin(),p.end()));for(int u=0;u<k;u++){s[1+u]=en[u]-st[u];for(int v=0;v<k;v++)s[1+u]+=c[u][v];for(int v=u+1;v<k;v++)s[1+k+idx[u][v]]=c[u][v]+c[v][u];}gr[s].push_back(q);}
cout<<"k="<<k<<" one_groups="<<gr.size()<<" total="<<(1<<ne)<<endl;vector<vector<int>>col;int maxsize=0;for(auto &kv:gr)if(kv.second.size()>1){col.push_back(kv.second);maxsize=max(maxsize,(int)kv.second.size());}cout<<"collision_groups="<<col.size()<<" max_size="<<maxsize<<endl;
for(int a=0;a<(1<<k)&&col.size();a++)for(int b=0;b<(1<<k)&&col.size();b++){vector<vector<int>>nc;for(auto &g:col){map<long long,vector<int>>gg;for(int q:g)gg[resp(outs(k,q),a,b)].push_back(q);for(auto &kv:gg)if(kv.second.size()>1)nc.push_back(kv.second);}if(nc.size()!=col.size())cout<<"probe="<<a<<","<<b<<" remaining_groups="<<nc.size()<<endl;col=nc;}
for(auto &g:col){cout<<"COLLISION";for(int x:g)cout<<" "<<x;cout<<endl;}cout<<"two_collisions="<<col.size()<<endl;
}
