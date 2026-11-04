#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using U=uint64_t;
vector<U> endpoints(const vector<int>&r){int n=r.size(),N=1<<n;vector<U>d(N*n);for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;for(int m=1;m<N;m++)for(int v=0;v<n;v++)if(m>>v&1){int rest=m^(1<<v);if(!rest)continue;U z=0;for(int u=0;u<n;u++)if((rest>>u&1)&&(r[u]>>v&1))z+=d[rest*n+u];d[m*n+v]=z;}return vector<U>(d.begin()+(N-1)*n,d.end());}
int module(const vector<int>&r){for(int m=1;m<511;m++){if(__builtin_popcount((unsigned)m)<2)continue;bool good=true;for(int v=0;v<9;v++)if(!(m>>v&1)&&(r[v]&m)!=0&&(r[v]&m)!=m){good=false;break;}if(good)return m;}return 0;}
int main(){string s;getline(cin,s);int node;U code;cout<<"node\tcode\tproper_module\tp_mask\tstarts\tends\n";while(cin>>node>>code){vector<int>r(9),op(9);int b=35;for(int j=1;j<9;j++)for(int i=0;i<j;i++,b--){if(code>>b&1)r[i]|=1<<j;else r[j]|=1<<i;}int m=module(r);cout<<node<<'\t'<<code<<'\t'<<m;if(m){cout<<"\t-\t-\t-\n";continue;}for(int i=0;i<9;i++)op[i]=511^r[i]^(1<<i);auto e=endpoints(r),a=endpoints(op);int p=0;for(int i=0;i<9;i++)p|=((a[i]+e[i])%2)<<i;cout<<'\t'<<p<<'\t';for(int i=0;i<9;i++){if(i)cout<<',';cout<<a[i];}cout<<'\t';for(int i=0;i<9;i++){if(i)cout<<',';cout<<e[i];}cout<<'\n';}}
