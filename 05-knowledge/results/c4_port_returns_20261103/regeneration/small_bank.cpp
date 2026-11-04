#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
using namespace std;using I=int64_t;using T=vector<unsigned>;
struct Profile{I h,p;vector<I>f,g,e,s;vector<vector<I>>j;};
Profile profile(const T&t){int n=t.size(),N=1<<n,full=N-1;vector<I>e(N*n),s(N*n),h(N);h[0]=1;
for(int mask=1;mask<N;mask++)for(unsigned bits=mask;bits;bits&=bits-1){int v=countr_zero(bits),rest=mask^(1<<v);I ev=rest?0:1,sv=ev;for(unsigned q=rest;q;q&=q-1){int u=countr_zero(q);if(t[u]>>v&1)ev+=e[rest*n+u];if(t[v]>>u&1)sv+=s[rest*n+u];}e[mask*n+v]=ev;s[mask*n+v]=sv;h[mask]+=ev;}
Profile out{h[full],0,vector<I>(n),vector<I>(n),vector<I>(n),vector<I>(n),vector<vector<I>>(n,vector<I>(n))};
for(int mask=0;mask<N;mask++){int other=full^mask;out.p+=h[mask]*h[other];for(int u=0;u<n;u++){out.f[u]+=e[mask*n+u]*h[other];out.g[u]+=s[mask*n+u]*h[other];for(int v=0;v<n;v++)out.j[u][v]+=e[mask*n+u]*s[other*n+v];}}
for(int v=0;v<n;v++){out.e[v]=e[full*n+v];out.s[v]=s[full*n+v];}return out;}
unsigned terminal(const T&t){int n=t.size();vector<unsigned>reach=t;unsigned all=(1u<<n)-1;for(int i=0;i<n;i++)reach[i]|=1u<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];for(auto row:reach)all&=row;return all;}
T birth(T t,unsigned cut){int n=t.size();for(int v=0;v<n;v++)if(cut>>v&1)t[v]|=1u<<n;t.push_back(((1u<<n)-1)^cut);return t;}
unsigned canonical(const T&t){int n=t.size();vector<int>order(n);for(int i=0;i<n;i++)order[i]=i;sort(order.begin(),order.end(),[&](int a,int b){return popcount(t[a])<popcount(t[b])||(popcount(t[a])==popcount(t[b])&&a<b);});vector<pair<int,int>>groups;for(int i=0;i<n;){int j=i+1;while(j<n&&popcount(t[order[i]])==popcount(t[order[j]]))j++;groups.push_back({i,j});i=j;}unsigned best=~0u;function<void(int)>go=[&](int g){if(g==(int)groups.size()){unsigned code=0;int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)if(t[order[i]]>>order[j]&1)code|=1u<<k;best=min(best,code);return;}auto[a,b]=groups[g];sort(order.begin()+a,order.begin()+b);do{go(g+1);}while(next_permutation(order.begin()+a,order.begin()+b));};go(0);return best;}
void print_hit(const T&t,const Profile&p,unsigned W,int c,int source,int cut){I f=0,g=0,j=0,e=0,s=0;for(int u=0;u<(int)t.size();u++)if(W>>u&1){f+=p.f[u];g+=p.g[u];e+=p.e[u];s+=p.s[u];for(int v=0;v<(int)t.size();v++)if(W>>v&1)j+=p.j[u][v];}
cout<<"{\"source\":"<<source<<",\"cut\":"<<cut<<",\"intercept\":"<<c<<",\"ports_mask\":"<<W<<",\"signature\":["<<p.h<<','<<p.p<<','<<f<<','<<g<<','<<j<<','<<e<<','<<s<<"],\"rows\":[";for(int i=0;i<(int)t.size();i++){if(i)cout<<',';cout<<t[i];}cout<<"]}";}
int main(int argc,char**argv){string mode=argc>1?argv[1]:"local";array<long long,2>counts{};int printed=0,base_representatives=0;vector<array<long long,3>>by_order;cout<<"{\"mode\":\""<<mode<<"\",\"hits\":[";
auto test=[&](const T&t,int source,int cut,bool keep_all){auto p=profile(t);auto term=terminal(t);array<int,2>goals=source==2?array<int,2>{-37,-35}:source==4?array<int,2>{-31,-29}:array<int,2>{-35,-29};
for(unsigned W=1;W<(1u<<t.size());W++){if(!(W&term))continue;I f=0,j=0;for(int u=0;u<(int)t.size();u++)if(W>>u&1){f+=p.f[u];for(int v=0;v<(int)t.size();v++)if(W>>v&1)j+=p.j[u][v];}if(f!=63)continue;for(int k=0;k<2;k++)if(p.h-j==goals[k]){bool first=counts[k]++==0;if(keep_all||first){if(printed++)cout<<',';print_hit(t,p,W,goals[k],source,cut);}}}};
if(mode=="complete7"){int n=6;vector<pair<int,int>>edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.push_back({i,j});map<unsigned,T>bank;for(unsigned code=0;code<(1u<<edges.size());code++){T t(n);for(int k=0;k<(int)edges.size();k++){auto[u,v]=edges[k];if(code>>k&1)t[u]|=1u<<v;else t[v]|=1u<<u;}bank.emplace(canonical(t),t);}base_representatives=bank.size();for(auto&[code,t]:bank)for(unsigned cut=0;cut<64;cut++)test(birth(t,cut),0,cut,true);by_order.push_back({7,counts[0],counts[1]});}
else if(mode=="local"||mode=="two"){for(int src:{2,4}){T t=src==2?T{30,28,8,16,4}:T{30,28,24,16,0};counts={0,0};for(unsigned cut=0;cut<32;cut++){auto X=birth(t,cut);if(mode=="local")test(X,src,cut,true);else for(unsigned second=0;second<64;second++)test(birth(X,second),src,cut+32*second,true);}by_order.push_back({src,counts[0],counts[1]});}}
else{for(int n=1;n<=6;n++){counts={0,0};vector<pair<int,int>>edges;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.push_back({i,j});for(unsigned code=0;code<(1u<<edges.size());code++){T t(n);for(int k=0;k<(int)edges.size();k++){auto[u,v]=edges[k];if(code>>k&1)t[u]|=1u<<v;else t[v]|=1u<<u;}test(t,0,-1,false);}by_order.push_back({n,counts[0],counts[1]});}}
cout<<"],\"counts\":[";for(int i=0;i<(int)by_order.size();i++){if(i)cout<<',';auto a=by_order[i];cout<<'['<<a[0]<<','<<a[1]<<','<<a[2]<<']';}cout<<"],\"base_representatives\":"<<base_representatives<<"}\n";return 0;}
