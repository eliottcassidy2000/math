#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
using V=vector<unsigned>; using Num=long long;
struct P{Num e=0,o=0; P operator+(P a)const{return {e+a.e,o+a.o};}
 P operator*(P a)const{return {e*a.e+o*a.o,e*a.o+o*a.e};}
 P operator*(Num a)const{return {e*a,o*a};}};
struct Tile{V r; Num h,R; vector<P> start,end; vector<vector<P>> through;};
Tile tile(const V&r){
 int n=r.size(),F=(1<<n)-1;
 vector<vector<vector<int>>> d(n,vector<vector<int>>(F+1,vector<int>(n)));
 vector<int> H(F+1);H[0]=1;
 for(int i=0;i<n;i++)d[i][1<<i][i]=1;
 for(int S=1;S<=F;S++)for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(d[i][S][j]){
  H[S]+=d[i][S][j];
  for(unsigned z=r[j]&~S;z;z&=z-1)d[i][S|(z&-z)][countr_zero(z)]+=d[i][S][j];
 }
 Tile t{r,H[F],0,vector<P>(n),vector<P>(n),vector<vector<P>>(n,vector<P>(n))};
 for(int S=0;S<=F;S++){
  t.R+=H[S]*H[F^S];
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   Num q=d[i][S][j]*H[F^S];P x=(popcount(unsigned(S))&1)?P{0,q}:P{q,0};
   t.start[i]=t.start[i]+x;t.end[j]=t.end[j]+x;t.through[i][j]=t.through[i][j]+x;
  }
 }
 return t;
}
unsigned canonical(const V&r){
 int n=r.size();V p(n);for(int i=0;i<n;i++)p[i]=i;
 unsigned z=~0u;
 do{unsigned c=0;for(int i=1;i<n;i++)for(int j=0;j<i;j++)c=2*c+(r[p[j]]>>p[i]&1);z=min(z,c);}
 while(next_permutation(p.begin(),p.end()));return z;
}
bool strong(const V&r){unsigned F=(1u<<r.size())-1;
 for(unsigned k=0;k<r.size();k++){
  unsigned b=0,S=1u<<k;while(b!=S){b=S;for(unsigned z=b;z;z&=z-1)S|=r[countr_zero(z)];}
  if(S!=F)return false;
 }return true;
}

map<array<int,4>,V> bank;
void search(V a){
 Tile t=tile(a);int n=a.size();
 if(t.h>41)return;
 if(strong(a))for(int i=0;i<n;i++){
    V b;for(int j=0;j<n;j++)if(i!=j)b.push_back((a[j]&((1u<<i)-1))|((a[j]>>(i+1))<<i));
    Num g=tile(b).h,h=t.h,E=t.start[i].e,O=t.start[i].o;
    if(O-E==g&&(h*h+(E+O)*(E+O)-g*g==1729)) {
      if(bank.emplace(array<int,4>{int(h),int(E),int(O),int(g)},a).second){
       cout<<"source "<<h<<" "<<E<<" "<<O<<" "<<g<<" root "<<i;
       for(auto x:a)cout<<" "<<x;cout<<"\n";
      }
    }
 }
 if(n==7)return;for(unsigned c=0;c<(1u<<(n-1));c++){
   V r=a;for(int i=0;i<n;i++)if(!(c>>i&1))r[i]|=1u<<n;
   r.push_back(c);search(r);
 }
}
int main(){search({0});cout<<"profiles "<<bank.size()<<"\n";if(!bank.empty())return 1;}
