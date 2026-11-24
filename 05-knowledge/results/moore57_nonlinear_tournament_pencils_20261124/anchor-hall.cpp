// Local Hall test for one or more regular-spread anchors in PG(3,7).
// Compile: c++ -O3 -std=c++17 anchor-hall.cpp -o anchor-hall
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <set>
#include <vector>
using namespace std;
constexpr int q=7, N=2401;
using M=array<int,4>;
using Mask=bitset<N>;
M V[N];
int inv[N];
int mod(int x){x%=q; return x<0?x+q:x;}
int id(M a){return ((a[0]*q+a[1])*q+a[2])*q+a[3];}
int mul(int a,int b){
 const auto [w,x,y,z]=V[a]; const auto [e,f,g,h]=V[b];
 return id({mod(w*e+x*g),mod(w*f+x*h),mod(y*e+z*g),mod(y*f+z*h)});
}
int add(int a,int b){
 const auto A=V[a],B=V[b];
 return id({mod(A[0]+B[0]),mod(A[1]+B[1]),mod(A[2]+B[2]),mod(A[3]+B[3])});
}
int sub(int a,int b){
 const auto A=V[a],B=V[b];
 return id({mod(A[0]-B[0]),mod(A[1]-B[1]),mod(A[2]-B[2]),mod(A[3]-B[3])});
}
vector<int> elems(const Mask &x){
 vector<int> v;for(int j=0;j<N;j++)if(x[j])v.push_back(j);return v;
}
int matching(const vector<vector<int>>& nbr,int n) {
 vector<int> L(n,-1),R(n,-1),d(n); int found=0;
 for(;;){
  queue<int> Q; fill(d.begin(),d.end(),-1);
  for(int i=0;i<n;i++)if(L[i]<0){d[i]=0;Q.push(i);}
  while(!Q.empty()){int i=Q.front();Q.pop();for(int j:nbr[i])if(R[j]>=0&&d[R[j]]<0){d[R[j]]=d[i]+1;Q.push(R[j]);}}
  function<bool(int)> dfs=[&](int u){
   for(int v:nbr[u])if(R[v]<0||(d[R[v]]==d[u]+1&&dfs(R[v]))){L[u]=v;R[v]=u;return true;}
   d[u]=-1;return false;
  };
  int add=0;for(int u=0;u<n;u++)if(L[u]<0&&dfs(u))add++;
  if(!add)return found;
  found+=add;
 }
}
int main(int argc,char**argv){
 vector<int> gl;
 for(int a=0;a<q;a++)for(int b=0;b<q;b++)for(int c=0;c<q;c++)for(int d=0;d<q;d++){
  int i=id({a,b,c,d});V[i]={a,b,c,d};inv[i]=-1;
  int det=mod(a*d-b*c);
  if(det){int z=1;while(mod(z*det)!=1)z++; inv[i]=id({mod(z*d),mod(-z*b),mod(-z*c),mod(z*a)});gl.push_back(i);}
 }
 set<vector<int>> fields;
 for(int i:gl){
  auto [a,b,c,d]=V[i];int det=mod(a*d-b*c);
  bool irreducible=true;for(int x=0;x<q;x++)if(!mod(x*x-(a+d)*x+det))irreducible=false;
  if(!irreducible)continue;
  vector<int> field;
  for(int s=0;s<q;s++)for(int t=0;t<q;t++)field.push_back(id({mod(s+t*a),mod(t*b),mod(t*c),mod(s+t*d)}));
  sort(field.begin(),field.end());fields.insert(field);
 }
 vector<Mask> spreads;
 for(const auto &F:fields){
  Mask seen;
  for(int a:gl)if(!seen[a]){
   Mask S;for(int f:F)S.set(mul(a,f));
   seen|=S;spreads.push_back(S);
  }
 }
 assert(fields.size()==21 && spreads.size()==882);
 auto S=spreads[0];
 vector<int> allowed;
 for(int i=0;i<int(spreads.size());i++)if((S&spreads[i]).count()==1)allowed.push_back(i);
 assert(allowed.size()==721);
 if(argc>=2 && string(argv[1])=="plane-slices"){
  int num_pass=0,tot_groups=0,min_bad=36;
  auto A=elems(S);
  auto lines=[&](const vector<int>& W){
   set<vector<int>> sets;
   for(int w:W)if(w){
    vector<int> vs;int x=w;
    for(int i=1;i<q;i++){vs.push_back(x);x=add(x,w);}
    sort(vs.begin(),vs.end());sets.insert(vs);
   }
   assert(sets.size()==8);
   return sets;
  };
  auto dl=lines(A);
  for(int j:allowed){
   Mask T=spreads[j];
   auto B=elems(T);
   // Inversion of nonzero matrices takes each scalar line to another
   // scalar line, so the right pinned slice ranges over lines in T.
   auto el=lines(B);
   bool bad[N][49]={}; int BID[N]={};
   for(int b=1;b<49;b++)BID[B[b]]=b;
   for(int D:A)if(D)for(int b=1;b<49;b++){
    int E=B[b],hits=0;
    for(int K:B){int x=add(D,K);if(inv[x]>=0&&S[sub(inv[x],inv[E])])hits++;}
    bad[D][b]=hits!=int(inv[sub(D,E)]>=0);
   }
   for(auto& U:dl)for(auto& V:el){
    int n=0;for(int D:U)for(int E:V)n+=bad[D][BID[E]];
    num_pass+=(n==0);min_bad=min(min_bad,n);tot_groups++;
   }
  }
  cout<<"adjacent plane-slice centers: "<<num_pass<<"/"<<tot_groups
      <<" pairs of forced zero slices pass; minimum failing pairs "<<min_bad<<endl;
  assert(num_pass==0 && tot_groups==721*64 && min_bad==12);
  return 0;
 }
 if(argc>=2 && string(argv[1])=="parallel"){
  auto A=elems(S);int ok=0, tested=0, minbad=48*48;
  for(int j:allowed){
   auto B=elems(spreads[j]);int bad=0;
   for(int D:A)if(D)for(int E:B)if(E){
    int hits=0;
    for(int K:B){
     int x=add(D,K);
     if(inv[x]>=0 && S[sub(inv[x],inv[E])])hits++;
    }
    if(hits!=int(inv[sub(D,E)]>=0))bad++;
   }
   minbad=min(minbad,bad);
   if(!bad)ok++;
   tested++;
  }
  assert(tested==721 && ok==0 && minbad==1752);
  cout<<"double-parallel anchors: "<<ok<<"/"<<tested<<"; minimum failing cross pairs "<<minbad<<endl;
  return 0;
 }
 int first=argc>=2?atoi(argv[1]):0,last=argc>=3?atoi(argv[2]):first+1;
 assert(first>=0&&last<=int(allowed.size())&&first<last);
 for(int k=first;k<last;k++){
  Mask T=spreads[allowed[k]],U=S|T;
  auto A=elems(S),B=elems(T);
  vector<int> X;
  for(int i:gl)if(!U[i])X.push_back(i);
  assert(X.size()==1920);
  vector<vector<int>> nbr;
  for(int D:A)if(D)for(int E:B)if(E && inv[sub(D,E)]>=0){
   // Transport D,E to horizontal,vertical. All finite forbidden anchor
   // lines skew to E can be tested as finite matrices M in this chart.
   auto transform=[&](int z){int ie=inv[sub(z,E)];return ie<0?-1:mul(sub(z,D),ie);};
   Mask forbid;for(int z:elems(U))if(z!=D && z!=E){int x=transform(z);if(x>=0)forbid.set(x);}
   forbid.set(id({1,0,0,1})); // original vertical line, not D or E
   Mask pass;
   for(auto &P:spreads)if((forbid&P).none())pass|=P;
   vector<int> ns;
   for(int i=0;i<int(X.size());i++){
    int x=transform(X[i]);if(x>=0 && x!=0 && pass[x])ns.push_back(i);
   }
   nbr.push_back(move(ns));
  }
  assert(nbr.size()==1920);
  int lo=1920,hi=0;long long edges=0;
  for(auto &ns:nbr){lo=min(lo,int(ns.size()));hi=max(hi,int(ns.size()));edges+=ns.size();}
  int m=matching(nbr,1920);
  assert(m==1920 && edges>=1080192 && edges<=1134624);
  cout<<"anchor "<<k<<": degree "<<lo<<".."<<hi<<", "<<edges<<" compatible incidences, Hall matching "<<m<<"/1920"<<endl;
 }
}
