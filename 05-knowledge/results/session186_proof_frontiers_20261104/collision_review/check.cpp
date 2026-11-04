#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;
using U=uint64_t;
U count(const vector<int>&r){int n=r.size();vector<int>in(n);for(int v=0;v<n;v++)for(int u=0;u<n;u++)if(r[u]>>v&1)in[v]|=1<<u;vector<U>d((1<<n)*n);for(int s=1;s<(1<<n);s++){int vs=s;while(vs){int v=__builtin_ctz((unsigned)vs);vs&=vs-1;int t=s^(1<<v);if(!t){d[s*n+v]=1;continue;}int ps=in[v]&t;U x=0;while(ps){int u=__builtin_ctz((unsigned)ps);ps&=ps-1;x+=d[t*n+u];}d[s*n+v]=x;}}return accumulate(d.end()-n,d.end(),U(0));}
bool pairfree(const vector<int>&r){int n=r.size();for(int a=0;a<n;a++)for(int b=a+1;b<n;b++){bool homogeneous=true;for(int w=0;w<n;w++)if(w!=a&&w!=b&&((r[w]>>a&1)!=(r[w]>>b&1))){homogeneous=false;break;}if(homogeneous)return false;}return true;}
int main(){int seed[3][6]={{16,5,33,39,46,3},{16,5,49,39,42,3},{18,4,49,39,42,3}};for(int n=11;n<=13;n++)for(int j=0;j<3;j++){vector<int> base(n);for(int u=0;u<n;u++)for(int v=u+1;v<n;v++){int su=min(u,5),sv=min(v,5);bool f=su==sv||((seed[j][su]>>sv)&1);base[f?u:v]|=1<<(f?v:u);}U h=count(base);vector<pair<int,int>>edges;for(int u=0;u<n;u++)for(int v=u+1;v<n;v++)edges.push_back({u,v});long total[4]={},free[4]={},same[4]={};int m=edges.size();for(int a=0;a<m;a++)for(int b=a;b<m;b++)for(int c=b;c<m;c++){int radius=(a==b?1:2)+(b==c?0:1);if(a==b&&b!=c)continue;vector<int>t=base;vector<int>ix={a};if(b!=a)ix.push_back(b);if(c!=b)ix.push_back(c);for(int e:ix){auto[u,v]=edges[e];t[u]^=1<<v;t[v]^=1<<u;}total[radius]++;if(!pairfree(t))continue;free[radius]++;if(count(t)==h){same[radius]++;cout<<"MATCH "<<n<<" "<<j;for(int x:t)cout<<" "<<x;cout<<"\n";}}for(int k=1;k<=3;k++)cout<<"ROW "<<n<<" "<<j<<" "<<k<<" "<<total[k]<<" "<<free[k]<<" "<<same[k]<<"\n";cout.flush();}}
