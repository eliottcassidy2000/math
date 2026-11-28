#include <cstdint>
#include <algorithm>
#include <array>
#include <cassert>
#include <bitset>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
constexpr int N=3250,D=56;
int n,k,delta;
vector<int> E[N];
bitset<N> A[N],P[N];
uint64_t choices;
U bit(int j){return U(1)<<j;}
int idx(U j){return __builtin_ctzll(j);}
void initialize(){
  for(int u=0;u<n;u++){
    bitset<N> B=A[u];
    for(int x:E[u]) B|=A[x];
    bitset<N> C=B;
    for(int x:E[u])for(int y:E[x]) C|=A[y];
    if((int)E[u].size()==k)continue;
    U br=0;
    for(int x:E[u])if(x>=k+1) br|=bit((x-k-1)/delta);
    for(int v=max(k+1,u+1);v<n;v++){
      if(u<k+1||(int)E[v].size()==k||C[v])continue;
      int i=(u-k-1)/delta,j=(v-k-1)/delta;
      if(i==j||(br&bit(j))) continue;
      bool ok=true;
      for(int x:E[v])if(x>=k+1&&(x-k-1)/delta==i){ok=false;break;}
      if(ok){P[u].set(v);P[v].set(u);choices++;}
    }
  }
}
bool augment(U* M,int* match,U &seen,int i){
  for(U z=M[i]&~seen;z;z&=z-1){
    int j=idx(z);
    if(seen&bit(j))continue;
    seen|=bit(j);
    if(match[j]<0||augment(M,match,seen,match[j])){match[j]=i;return true;}
  }
  return false;
}
int main(){
 cin>>n; assert(n<=N&&n>=5);
 for(int u=0;u<n;u++){
   int m;cin>>m;
   while(m--){int v;cin>>v;A[u].set(v);E[u].push_back(v);}
 }
 k=E[0].size();delta=k-1;
 assert(n==k*k+1&&delta<=D);
 uint64_t edge_pentagons=0;
 for(int u=0;u<n;u++)for(int v:E[u])if(u<v)
   for(int x:E[u])if(x!=v)for(int y:E[v])if(y!=u)
     edge_pentagons+=(A[x]&A[y]).count();
 assert(edge_pentagons%5==0);
 cout<<"Pentagons "<<edge_pentagons/5<<" vs completion "
     <<uint64_t(n)*k*delta*delta/10<<endl;
 initialize();
 cout<<"Initial optional edges "<<choices<<endl;
 vector<int> F;for(int u=0;u<n;u++)if((int)E[u].size()==k)F.push_back(u);
 for(int pass=1;;pass++){
   int boards=0,totalFixed=0,localForced=0,fail=0,weakest=D*D+1;
   uint64_t deleted=0, options=0;
   set<pair<int,int>> forced;
   for(int u:F)for(int v:F)if(u<v&&!A[u][v]){
      auto mid=A[u]&A[v];
      if(mid.count()!=1){cout<<"Missing full wedge "<<u<<","<<v<<endl;return 1;}
      int w=0;while(!mid[w])w++;
      vector<int> L,R;
      for(int x:E[u])if(x!=w)L.push_back(x);
      for(int y:E[v])if(y!=w)R.push_back(y);
      U M[D]{},S[D]{};
      int m=L.size(),fixed=0;
      assert(m==(int)R.size());
      for(int i=0;i<m;i++)for(int j=0;j<m;j++){
         if(A[L[i]][R[j]]){M[i]|=bit(j);fixed++;}
         if(P[L[i]][R[j]]) M[i]|=bit(j);
      }
      boards++;totalFixed+=fixed; int z=0;
      for(int i=0;i<m;i++)z+=__builtin_popcountll(M[i]);
      options+=z-fixed; if(fixed<m)weakest=min(weakest,z-fixed);
      int match[D];fill(match,match+D,-1);
      bool ok=true;
      for(int i=0;i<m;i++){U seen=0;if(!augment(M,match,seen,i)){ok=false;break;}}
      if(!ok){if(++fail<=10)cout<<"Wedge Hall failure "<<u<<","<<w<<","<<v<<endl;continue;}
      for(int i=0;i<m;i++)for(U Z=M[i];Z;Z&=Z-1)S[i]|=bit(match[idx(Z)]);
      for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(S[j]&bit(i)) S[j]|=S[i];
      for(int i=0;i<m;i++){
        int supported=0,last=-1;
        for(U Z=M[i];Z;Z&=Z-1){
           int j=idx(Z);
           if(match[j]==i||(S[match[j]]&bit(i))){supported++;last=j;continue;}
           assert(!A[L[i]][R[j]]);
           if(P[L[i]][R[j]]){P[L[i]].reset(R[j]);P[R[j]].reset(L[i]);deleted++;}
        }
        if(supported==1&&!A[L[i]][R[last]]){
           localForced++;
           forced.insert(minmax(L[i],R[last]));
        }
      }
   }
   cout<<"Pass "<<pass<<" boards "<<boards<<" fixed incidences "<<totalFixed
       <<" optional incidences "<<options<<" least on unfinished "<<(weakest>D*D?-1:weakest)
       <<" pentagon wedge deficit "<<uint64_t(boards)*delta-totalFixed
       <<" failures "<<fail<<" deleted "<<deleted<<" forced "<<forced.size()
       <<" (local incidences "<<localForced<<")"<<endl;
   if(fail)return 1;
   if(!deleted)break;
 }
}
