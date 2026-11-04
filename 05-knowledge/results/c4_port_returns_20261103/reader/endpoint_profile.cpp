// Fresh literal endpoint and ordered-joint profiles for the actual named graphs.
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;
using I=int64_t;
void need(bool ok,const char*m){if(!ok){cerr<<m<<endl;exit(2);}}
void print_array(const vector<I>&v){cout<<'[';for(int i=0;i<(int)v.size();i++){if(i)cout<<',';cout<<v[i];}cout<<']';}
int main(){
 int id,n;
 while(cin>>id>>n){
  need(n>=1&&n<=17,"declared finite profile order");
  I bound=1;for(int k=2;k<=n+1;k++){need(bound<=numeric_limits<I>::max()/k,"two-slot factorial bound");bound*=k;}
  vector<unsigned>t(n);for(auto&row:t)cin>>row;
  unsigned full=(1u<<n)-1;
  for(int i=0;i<n;i++){
   need(!(t[i]>>i&1)&&!(t[i]&~full),"valid row bits");
   for(int j=i+1;j<n;j++)need(((t[i]>>j&1)+(t[j]>>i&1))==1,"actual tournament");
  }
  int N=1<<n;vector<I>e(N*n),s(N*n),h(N);h[0]=1;
  for(int mask=1;mask<N;mask++)for(unsigned bits=mask;bits;bits&=bits-1){
   int v=countr_zero(bits),rest=mask^(1<<v);I ev=rest?0:1,sv=ev;
   for(unsigned q=rest;q;q&=q-1){int u=countr_zero(q);if(t[u]>>v&1)ev+=e[rest*n+u];if(t[v]>>u&1)sv+=s[rest*n+u];}
   e[mask*n+v]=ev;s[mask*n+v]=sv;h[mask]+=ev;
  }
  vector<I>f(n),g(n),ee(n),ss(n);vector<vector<I>>j(n,vector<I>(n));I p=0;
  for(int mask=0;mask<N;mask++){
   int other=full^mask;p+=h[mask]*h[other];
   for(int u=0;u<n;u++){
    I end=e[mask*n+u],start=s[mask*n+u];f[u]+=end*h[other];g[u]+=start*h[other];
    if(end)for(int v=0;v<n;v++)j[u][v]+=end*s[other*n+v];
   }
  }
  need(p<=bound,"P2 factorial bound");
  for(int v=0;v<n;v++){ee[v]=e[full*n+v];ss[v]=s[full*n+v];need(j[v][v]==0,"disjoint joint diagonal");}
  cout<<"{\"id\":"<<id<<",\"H\":"<<h[full]<<",\"P2\":"<<p<<",\"E1\":";print_array(ee);
  cout<<",\"S1\":";print_array(ss);cout<<",\"E2\":";print_array(f);cout<<",\"S2\":";print_array(g);cout<<",\"J\":[";
  for(int u=0;u<n;u++){if(u)cout<<',';print_array(j[u]);}cout<<"]}\n";
 }
}
