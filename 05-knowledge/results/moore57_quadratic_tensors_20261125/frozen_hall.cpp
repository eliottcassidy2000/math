// Branch-pair perfect matching support, followed by frozen radius-two clauses.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
using W=uint64_t;
constexpr int N=3250,K=57,D=56,F=58;
vector<int>A[N];
W opt[N][K];
int leaf(int b,int x){return F+b*D+x;}
int br(int u){return (u-F)/D;}
int sh(int u){return (u-F)%D;}
W bit(int x){return 1ULL<<x;}
void distance(int u,int *r){
 fill(r,r+N,4);r[u]=0;
 vector<int>Q{u};
 for(size_t i=0;i<Q.size();i++){
  int v=Q[i];if(r[v]==3)continue;
  for(int w:A[v])if(r[w]>r[v]+1){r[w]=r[v]+1;Q.push_back(w);}
 }
}
bool has(int u,int v){return v>=F&&((opt[u][br(v)]>>sh(v))&1);}
int match[D]; W seen;
bool augment(const W* X,int i){
 for(W S=X[i]&~seen;S;S&=S-1){
  int j=__builtin_ctzll(S);
  if(seen&bit(j))continue;
  seen|=bit(j);
  if(match[j]<0||augment(X,match[j])){match[j]=i;return true;}
 }
 return false;
}
int main(){
 int n;cin>>n;if(n!=N)return 2;
 for(auto& X:A){int k;cin>>k;while(k--){int u;cin>>u;X.push_back(u);}}
 vector<int>full;
 long long cnt=0,removed=0;
 int R[N];
 for(int u=0;u<N;u++)if(A[u].size()==57)full.push_back(u);
 // Directed copy of eligibility. A fixed branch slot at either endpoint
 // disallows all other contacts in that branch.
 W occupied[N]{};
 for(int u=F;u<N;u++)for(int w:A[u])if(w>=F)occupied[u]|=bit(br(w));
 for(int u=F;u<N;u++)if(A[u].size()<57){
  distance(u,R);
  for(int v=u+1;v<N;v++)if(br(u)!=br(v)&&A[v].size()<57&&
     !(occupied[u]&bit(br(v)))&&!(occupied[v]&bit(br(u)))&&R[v]>3){
   opt[u][br(v)]|=bit(sh(v));opt[v][br(u)]|=bit(sh(u));++cnt;
  }
 }
 cout<<"full vertices "<<full.size()<<"; eligible undirected edges "<<cnt<<endl;
 int hallFail=0;
 for(int i=0;i<K;i++)for(int j=i+1;j<K;j++){
  W X[D]{},free=0;
  for(int y=0;y<D;y++){
   int u=leaf(i,y);
   X[y]=opt[u][j];
   if(!(occupied[u]&bit(j)))free|=bit(y);
  }
  if(!free)continue;
  fill(match,match+D,-1);
  bool ok=true;
  for(W U=free;U;U&=U-1){
   seen=0;
   if(!augment(X,__builtin_ctzll(U))){ok=false;break;}
  }
  if(!ok){if(++hallFail<=8)cout<<"Matching failure "<<i<<","<<j<<endl;continue;}
  W reach[D]{};
  for(int y=0;y<D;y++)for(W Z=X[y];Z;Z&=Z-1)
   { int u=match[__builtin_ctzll(Z)]; assert(u>=0); reach[y]|=bit(u); }
  for(int k=0;k<D;k++)for(int y=0;y<D;y++)if(reach[y]&bit(k))reach[y]|=reach[k];
  for(W U=free;U;U&=U-1){
   int y=__builtin_ctzll(U),u=leaf(i,y);
   for(W Z=X[y];Z;Z&=Z-1){
    int z=__builtin_ctzll(Z);
    if(match[z]!=y&&!(reach[match[z]]&bit(y))){
     opt[u][j]&=~bit(z);opt[leaf(j,z)][i]&=~bit(y);++removed;
    }
   }
  }
 }
 cout<<"Hall failures "<<hallFail<<"; nonextendible matching edges "<<removed<<endl;
 std::set<pair<int,int>> units;
 long long failures=0,checks=0;
 array<long long,58> hist{};
 for(int v:full){
  distance(v,R);
  for(int u=F;u<N;u++)if(R[u]>2){
   ++checks;int last=-1,num=0;
   for(int w:A[v])if(has(u,w)){last=w;++num;}
   ++hist[num];
   if(!num){if(++failures<=8)cout<<"Radius demand empty "<<v<<","<<u<<endl;}
   if(num==1)units.insert(minmax(u,last));
  }
 }
 cout<<"frozen radius demands "<<checks<<"; failures "<<failures<<"; forced edges "<<units.size()<<endl;
 cout<<"signature";for(int k=0;k<58;k++)if(hist[k])cout<<" "<<k<<":"<<hist[k];cout<<endl;
 for(auto [u,v]:units)cout<<"Force "<<u<<" "<<v<<endl;
}
