#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
int n,nn,k; vector<pair<int,int>> F;
vector<int> owner; vector<int> H;
array<int,9>A{},U{};
void paths(int i,int S,int c){
 if(i==n){
  unsigned t=0;
  for(int z=0;z<k;z++)if((A[U[F[z].first]]>>U[F[z].second])&1)t|=1u<<z;
  if(owner[t]!=-1&&owner[t]!=c)throw 10;
  owner[t]=c; ++H[c];return;
 }
 unsigned d=(i?A[U[i-1]]:nn)&~S;
 while(d){unsigned b=d&-d;d^=b;U[i]=__builtin_ctz(b);paths(i+1,S|b,c);}
}
vector<int> dsu;int find(int i){return dsu[i]==i?i:dsu[i]=find(dsu[i]);}
void join(int i,int j){i=find(i);j=find(j);if(i>j)swap(i,j);dsu[j]=i;}

int main(int argc,char**argv){
 if(argc!=2)return 2;
 n=stoi(argv[1]);if(n<3||n>8)return 2;
 nn=(1<<n)-1;
 for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)F.emplace_back(i,j);
 k=F.size();unsigned cap=1u<<k;
 owner.assign(cap,-1);
 string s;int C=0;
 while(getline(cin,s)){
  if(s.empty()||s[0]<'0'||s[0]>'9')continue;
  istringstream st(s);
  for(int i=0;i<n;i++)if(!(st>>A[i]))return 2;
  H.push_back(0); paths(0,0,C);++C;
 }
 vector<unsigned> mirror(k);
 for(int j=0;j<k;j++){
  auto p=make_pair(n-1-F[j].second,n-1-F[j].first);
  mirror[j]=1u<<(find(F.begin(),F.end(),p)-F.begin());
 }
 auto reflect=[&](unsigned j){unsigned r=0;while(j){unsigned b=j&-j;j^=b;r|=mirror[__builtin_ctz(b)];}return r;};
 dsu.resize(C);iota(dsu.begin(),dsu.end(),0);
 for(unsigned t=0;t<cap;t++)if(owner[t]<0)return 4;
 for(unsigned t=0;t<cap;t++)join(owner[t],owner[reflect(t)]);
 vector<uint64_t> E;vector<int> rootId(C,-1),val;
 for(int i=0;i<C;i++)if(find(i)==i){rootId[i]=val.size();val.push_back(H[i]);}
 vector<int> ids(cap);
 for(unsigned t=0;t<cap;t++){
  int r=find(owner[t]);if(H[r]!=H[owner[t]])return 4;
  ids[t]=rootId[r];
 }
 auto add=[&](int a,int b,bool blue){
  if(a>b)swap(a,b);
  E.push_back(((uint64_t)a*val.size()+b)*2+blue);
 };
 for(unsigned t=0;t<cap/2;t++)add(ids[t],ids[cap-1-t],reflect(t)==t);
 sort(E.begin(),E.end());E.erase(unique(E.begin(),E.end()),E.end());
 vector<vector<int>> G(val.size());int edges=0,loops=0,BE=0;
 for(unsigned j=0;j<E.size();){
  uint64_t e=E[j]/2;int a=e/val.size(),b=e%val.size();
  bool blue=(E[j]&1)!=0; ++j;
  if(j<E.size()&&E[j]/2==e){blue=true;++j;}
  if(a==b){loops++;continue;}
  edges++;BE+=blue;
  G[a].push_back(b);G[b].push_back(a);
 }
 vector<int>D(val.size(),-1);int r=ids.back();queue<int>Q;Q.push(r);D[r]=0;
 while(!Q.empty()){int a=Q.front();Q.pop();for(int b:G[a])if(D[b]<0){D[b]=D[a]+1;Q.push(b);}}
 int K=*max_element(D.begin(),D.end());vector<int>S(K+1);vector<set<int>>hs(K+1);
 for(unsigned v=0;v<D.size();v++)if(D[v]>=0){S[D[v]]++;hs[D[v]].insert(val[v]);}
 cout<<"n "<<n<<" classes "<<C<<" vertices "<<val.size()<<" edges "<<edges<<" blue-support "<<BE<<" loops "<<loops<<" reachable "<<accumulate(S.begin(),S.end(),0)<<"\n";
 cout<<"shell vertices";for(int c:S)cout<<' '<<c;cout<<"\nshell distinct H";for(auto&t:hs)cout<<' '<<t.size();cout<<'\n';
}
