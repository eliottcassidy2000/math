#include <bits/stdc++.h>
using namespace std;using R=vector<unsigned>;using Z=long long;
void need(bool b,string s){if(!b)throw runtime_error(s);}
Z pw(Z x,int n){Z r=1;while(n--)r*=x;return r;}
// Independent labels: 0=a, 1=u, 2=b; proof order, not source implementation order.
R make_graph(int i,int k,int j,int m){vector<int>w{1};w.insert(w.end(),i,0);w.push_back(1);w.insert(w.end(),j,0);w.push_back(2);w.insert(w.end(),k,0);w.push_back(2);int q=w.size(),n=3+q+m;R r(n);auto arc=[&](int a,int b){r[a]|=1u<<b;};for(int a=0;a<3;a++)for(int b=a+1;b<3;b++)arc(a,b);for(int t=0;t<q;t++){int v=3+t;for(int x=0;x<3;x++){bool oldbeats=(x==0&&w[t]==1)||(x==2&&w[t]==2);if(oldbeats)arc(x,v);else arc(v,x);}for(int t2=t+1;t2<q;t2++)arc(v,3+t2);for(int x=3+q;x<n;x++)arc(v,x);}for(int x=3+q;x<n;x++){for(int y=0;y<3;y++)arc(x,y);for(int y=x+1;y<n;y++)arc(x,y);}return r;}
// Inclusion-exclusion over allowed vertices of walks of exactly n-1 edges.
// No Hamilton endpoint-subset recurrence or Hamilton-word enumeration is used.
Z walk_ie(R const&r){int n=r.size(),N=1<<n;Z ans=0;vector<Z>d(n),e(n);for(int S=1;S<N;S++){fill(d.begin(),d.end(),0);for(int a=0;a<n;a++)if(S>>a&1)d[a]=1;for(int t=1;t<n;t++){fill(e.begin(),e.end(),0);for(int a=0;a<n;a++)if(d[a]){unsigned bits=r[a]&S;while(bits){int b=__builtin_ctz(bits);bits&=bits-1;e[b]+=d[a];}}d.swap(e);}Z h=accumulate(d.begin(),d.end(),0LL);ans+=((n-__builtin_popcount((unsigned)S))%2?-h:h);}return ans;}
bool strong(R const&r){int n=r.size();for(int x=0;x<n;x++){unsigned z=1u<<x,old=0;while(old!=z){old=z;for(int a=0;a<n;a++)if(z>>a&1)z|=r[a];}if(z!=(1u<<n)-1)return false;}return true;}
R birth(R r,bool left,vector<int>const&exception){int n=r.size();r.push_back(0);for(int a=0;a<n;a++){bool ex=find(exception.begin(),exception.end(),a)!=exception.end();bool newbeats=left?!ex:ex;if(newbeats)r[n]|=1u<<a;else r[a]|=1u<<n;}return r;}
void check_map(R const&a,R const&b,vector<int>const&f){need(a.size()==b.size(),"map order");for(int x=0;x<(int)a.size();x++)for(int y=0;y<(int)a.size();y++)need(bool(a[x]>>y&1)==bool(b[f[x]]>>f[y]&1),"literal physical mapping");}
int main(){int grids=0,returns=0;vector<array<Z,5>>counts;for(int i=0;i<=1;i++)for(int k=0;k<=1;k++)for(int j=0;j<=2;j++)for(int m=0;m<=3;m++){R r=make_graph(i,k,j,m);Z expected=4*(3+2*i)*(5+8*k)*pw(2,j)*pw(3,m)+(j?1:-1);need(walk_ie(r)==expected,"walk tensor count");if(m)need(strong(r),"strong grid");grids++;counts.push_back({i,k,j,m,expected});}
// Every word of length <=5 on the two returning moves, preserving the same old ports.
for(int i=0;i<=1;i++)for(int k=0;k<=1;k++)for(int len=0;len<=5;len++)for(int code=0;code<(1<<len);code++){
 int j=1,m=1;R T=make_graph(i,k,j,m);vector<int>oldToNew(T.size());iota(oldToNew.begin(),oldToNew.end(),0);vector<int>leftports;for(int t=0;t<2+i;t++)leftports.push_back(3+t);vector<int>middle{5+i};vector<int>roots{(int)T.size()-1};
 for(int step=0;step<len;step++){bool L=code>>step&1;int n=T.size();T=birth(T,L,L?leftports:vector<int>{0,1,2});if(L){middle.insert(middle.begin(),n);j++;}else{roots.push_back(n);m++;}}
 R U=make_graph(i,k,j,m);vector<int>f(T.size());f[0]=0;f[1]=1;f[2]=2;int next=3;for(int t=0;t<2+i;t++)f[3+t]=next++;for(int x:middle)f[x]=next++;int base=6+i;for(int t=0;t<2+k;t++)f[base+t]=next++;for(int x:roots)f[x]=next++;need(next==(int)T.size(),"map complete");check_map(T,U,f);returns++;
 }
// First insertion at the empty-middle boundary is the same left birth, not 2H-1.
for(int i=0;i<=1;i++)for(int k=0;k<=1;k++)for(int m=1;m<=3;m++){
 R T=make_graph(i,k,0,m);vector<int>E;for(int t=0;t<2+i;t++)E.push_back(3+t);R V=birth(T,true,E),U=make_graph(i,k,1,m);vector<int>f(V.size());int cut=5+i;for(int x=0;x<(int)T.size();x++)f[x]=x<cut?x:x+1;f.back()=cut;check_map(V,U,f);need(walk_ie(V)==2*walk_ie(T)+3,"empty boundary");returns++;
 }
cout<<"{\"status\":\"PASS\",\"walk_grid_graphs\":"<<grids<<",\"literal_port_word_and_boundary_maps\":"<<returns<<",\"grid\":[";for(int z=0;z<(int)counts.size();z++){if(z)cout<<',';cout<<'[';for(int x=0;x<5;x++){if(x)cout<<',';cout<<counts[z][x];}cout<<']';}cout<<"]}\n";}
