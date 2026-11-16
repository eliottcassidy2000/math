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
map<int,map<pair<int,unsigned>,Tile>> pool;
void collect(V r){
 Tile t=tile(r);if(t.h>15)return;
 if((t.h==3||t.h==5||t.h==15)&&strong(r))pool[t.h].emplace(make_pair(int(r.size()),canonical(r)),t);
 if(r.size()==6)return;
 unsigned n=r.size();
 for(unsigned z=0;z<(1u<<(n-1));z++){
  V a=r;for(unsigned i=0;i<n;i++)if(!(z>>i&1))a[i]|=1u<<n;
  a.push_back(z);collect(a);
 }
}
constexpr int BUF=32;
map<pair<Num,Num>,array<V,3>> menus;
Num N, target, stop_at=0;
int maxn,layouts_seen=0;
vector<pair<Num,Num>> triples;
void inspect(const vector<const Tile*>& a){
 int l=a.size(),n=0;int pos[BUF]={};Num pre[BUF]={1},suf[BUF]={};suf[l]=1;
 for(int i=0;i<l;i++){pre[i+1]=pre[i]*a[i]->h;n+=a[i]->r.size();pos[i+1]=n;}
 for(int i=l-1;i>=0;i--)suf[i]=suf[i+1]*a[i]->h;
 if(pre[l]!=15||n>=BUF)abort();
 layouts_seen++;
 Num W[BUF][BUF]={};P root[BUF]={};
 for(int i=0;i<l;i++){
  P tail{1,0};for(int j=i+1;j<l;j++)tail=tail*P{a[j]->R/2,a[j]->R/2};
  for(int u=pos[i];u<pos[i+1];u++){
   root[u]=a[i]->start[u-pos[i]]*tail*pre[i];
   for(int v=pos[i];v<pos[i+1];v++)W[u][v]=a[i]->through[u-pos[i]][v-pos[i]].e*pre[i]*suf[i+1];
  }
  P mid{1,0};
  for(int j=i+1;j<l;j++){
   for(int u=pos[i];u<pos[i+1];u++)for(int v=pos[j];v<pos[j+1];v++)
    W[u][v]=(a[i]->start[u-pos[i]]*mid*a[j]->end[v-pos[j]]).e*pre[i]*suf[j+1];
   mid=mid*P{a[j]->R/2,a[j]->R/2};
  }
 }
 const auto record=[&](unsigned A,Num d){
  if(d<target||d>target+2)return;
  unsigned F=(1u<<n)-1;
  if(!(A&((1u<<pos[1])-1))||!((F^A)&(F^((1u<<pos[l-1])-1))))return;
  P t{0,15};for(unsigned z=A;z;z&=z-1){P v=root[countr_zero(z)];t=t+P{v.o,v.e};}
  auto& ar=menus[{t.e,t.o}];
  if(!ar[d-target].empty())return;
  V r(n+1);
  for(int i=0;i<l;i++){
   unsigned later=(F^((1u<<pos[i+1])-1));
   for(int u=pos[i];u<pos[i+1];u++)
    r[u]=(a[i]->r[u-pos[i]]<<pos[i])|later|((!(A>>u&1))<<n);
  }
  r[n]=A;ar[d-target]=r;
  if(all_of(ar.begin(),ar.end(),[](const V&v){return !v.empty();})){
   triples.push_back({t.e,t.o});
   cout<<"TRIPLE "<<N<<" "<<t.e<<" "<<t.o<<" 15"<<endl;
   for(int i=0;i<3;i++){cout<<N+2*i;for(unsigned x:ar[i])cout<<" "<<x;cout<<endl;}
  }
 };
 // Left-to-right cut assignment; each edge is charged as soon as both bits are assigned.
 auto walk=[&](auto&&self,int i,unsigned A,unsigned B,Num d)->void{
  if(d>target+2 || int(triples.size())==stop_at)return;
  if(i==n){record(A,d);return;}
  Num x=d,y=d;
  for(unsigned t=A;t;t&=t-1)y+=W[countr_zero(t)][i];
  for(unsigned t=B;t;t&=t-1)x+=W[i][countr_zero(t)];
  self(self,i+1,A|(1u<<i),B,x);
  self(self,i+1,A,B|(1u<<i),y);
 };
 walk(walk,0,0,0,0);
}
void layout(const vector<const Tile*>& parts,int k,vector<const Tile*>& q,int depth,const Tile* one){
 if(stop_at==int(triples.size()))return;
 if(k==int(parts.size())){
  size_t p=q.size();
  for(int s=depth;s<=maxn;s++){inspect(q);q.push_back(one);if(stop_at==int(triples.size()))break;}
  q.resize(p);return;
 }
 auto t=parts[k];size_t p=q.size();
 for(int d=depth;d+int(t->r.size())<=maxn;d++){
  q.push_back(t);layout(parts,k+1,q,d+t->r.size(),one);q.pop_back();q.push_back(one);
  if(stop_at==int(triples.size()))break;
 }
 q.resize(p);
}
int main(int argc,char**argv){
 if(argc!=4)return 2;
 N=atoll(argv[1]);maxn=atoi(argv[2]);stop_at=atoi(argv[3]);target=(N-15)/2;
 if(N<15||N%2==0||maxn<5||maxn>22||stop_at<1)return 2;
 collect({0});
 if(pool[3].size()!=1||pool[5].size()!=1||pool[15].size()!=4)abort();
 auto c=&pool[3].begin()->second,b=&pool[5].begin()->second;Tile s=tile({0});
 vector<const Tile*> q;
 for(auto&[k,t]:pool[15])layout({&t},0,q,0,&s);
 layout({c,b},0,q,0,&s);layout({b,c},0,q,0,&s);
 cout<<"LAYOUTS "<<layouts_seen<<" PROFILES "<<menus.size()<<" TRIPLES "<<triples.size()<<endl;
}
