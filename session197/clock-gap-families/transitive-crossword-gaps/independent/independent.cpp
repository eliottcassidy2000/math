#include <bits/stdc++.h>
using namespace std;using R=vector<unsigned>;using Z=long long;
void need(bool b,string s){if(!b)throw runtime_error(s);}
Z pw(Z x,int n){Z r=1;while(n--)r*=x;return r;}
// Independent labels: 0=a, 1=u, 2=b; proof order, not source implementation order.
R make_graph(int i,int k,int j,int m){vector<int>w{1};w.insert(w.end(),i,0);w.push_back(1);w.insert(w.end(),j,0);w.push_back(2);w.insert(w.end(),k,0);w.push_back(2);int q=w.size(),n=3+q+m;R r(n);auto arc=[&](int a,int b){r[a]|=1u<<b;};for(int a=0;a<3;a++)for(int b=a+1;b<3;b++)arc(a,b);for(int t=0;t<q;t++){int v=3+t;for(int x=0;x<3;x++){bool oldbeats=(x==0&&w[t]==1)||(x==2&&w[t]==2);if(oldbeats)arc(x,v);else arc(v,x);}for(int t2=t+1;t2<q;t2++)arc(v,3+t2);for(int x=3+q;x<n;x++)arc(v,x);}for(int x=3+q;x<n;x++){for(int y=0;y<3;y++)arc(x,y);for(int y=x+1;y<n;y++)arc(x,y);}return r;}
// Inclusion-exclusion over allowed vertices of walks of exactly n-1 edges.
// No Hamilton endpoint-subset recurrence or Hamilton-word enumeration is used.
Z walk_ie(R const&r){int n=r.size(),N=1<<n;__int128 ans=0;vector<__int128>d(n),e(n);for(int S=1;S<N;S++){fill(d.begin(),d.end(),0);for(int a=0;a<n;a++)if(S>>a&1)d[a]=1;for(int t=1;t<n;t++){fill(e.begin(),e.end(),0);for(int a=0;a<n;a++)if(d[a]){unsigned bits=r[a]&S;while(bits){int b=__builtin_ctz(bits);bits&=bits-1;e[b]+=d[a];}}d.swap(e);}__int128 h=accumulate(d.begin(),d.end(),(__int128)0);ans+=((n-__builtin_popcount((unsigned)S))%2?-h:h);}need(ans>=0&&ans<LLONG_MAX,"walk result range");return (Z)ans;}
bool strong(R const&r){int n=r.size();for(int x=0;x<n;x++){unsigned z=1u<<x,old=0;while(old!=z){old=z;for(int a=0;a<n;a++)if(z>>a&1)z|=r[a];}if(z!=(1u<<n)-1)return false;}return true;}
R birth(R r,bool left,vector<int>const&exception){int n=r.size();r.push_back(0);for(int a=0;a<n;a++){bool ex=find(exception.begin(),exception.end(),a)!=exception.end();bool newbeats=left?!ex:ex;if(newbeats)r[n]|=1u<<a;else r[a]|=1u<<n;}return r;}
void check_map(R const&a,R const&b,vector<int>const&f){need(a.size()==b.size(),"map order");for(int x=0;x<(int)a.size();x++)for(int y=0;y<(int)a.size();y++)need(bool(a[x]>>y&1)==bool(b[f[x]]>>f[y]&1),"literal physical mapping");}
int main(){vector<array<Z,8>> out;for(int i=0;i<=2;i++)for(int k=0;k<=3;k++)for(int j=0;j<=2;j++)for(int m=0;m<=3;m++){
 Z Ui=i?pw(2,i)-1:1,Uk=k?pw(2,k)-1:1,f=pw(2,i+1)+1;
 Z A=4*f*(pw(3,k+1)+pw(2,k+1))*pw(2,j),D=-2*(f*(pw(3,k)-pw(2,k+1)+1)*pw(2,j)+(pw(2,k+1)+1)*(Ui-1)),C=Ui*Uk-2*(j==0);
 R T=make_graph(i,k,j,m);Z H=walk_ie(T);need(H==A*pw(3,m)+D*pw(2,m)+C,"independent 144-grid value");if(m)need(strong(T),"positive-root strongness");need((D==0)==(i<=1&&k<=1),"cancellation grid");out.push_back({i,k,j,m,A,D,C,H});}
 cout<<"{\"status\":\"PASS\",\"independent_walk_grid_graphs\":"<<out.size()<<",\"grid\":[";for(int q=0;q<(int)out.size();q++){if(q)cout<<',';cout<<'[';for(int z=0;z<8;z++){if(z)cout<<',';cout<<out[q][z];}cout<<']';}cout<<"]}\n";}
