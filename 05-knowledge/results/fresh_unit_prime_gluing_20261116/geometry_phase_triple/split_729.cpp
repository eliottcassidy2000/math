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

using Key=uint64_t;
Key canon(const V&r){
 int n=r.size();V p(n);
 array<V,12> groups;
 for(int i=0;i<n;i++)groups[popcount(r[i])].push_back(i);
 Key best=~Key(0);
 auto place=[&](auto&& self,int d,int at)->void{
  while(d<n&&groups[d].empty())d++;
  if(d==n){Key key=0;for(int j=1;j<n;j++)for(int i=0;i<j;i++)key=2*key+((r[p[i]]>>p[j])&1);best=min(best,key);return;}
  auto& a=groups[d];
  do{copy(a.begin(),a.end(),p.begin()+at);self(self,d+1,at+a.size());}while(next_permutation(a.begin(),a.end()));
 };
 place(place,0,0);return best;
}
bool strong(const V&r){unsigned F=(1u<<r.size())-1;
 for(unsigned k=0;k<r.size();k++){
  unsigned b=0,S=1u<<k;while(b!=S){b=S;for(unsigned z=b;z;z&=z-1)S|=r[countr_zero(z)];}
  if(S!=F)return false;
 }return true;
}
Num quick(const V&a){
 int n=a.size(),F=(1<<n)-1,d[2048][11];
 for(int S=1;S<=F;S++)for(unsigned z=S;z;z&=z-1){
  int i=countr_zero(z),T=S^(1<<i),h=!T;
  for(unsigned w=T&a[i];w;w&=w-1)h+=d[T][countr_zero(w)];
  d[S][i]=h;
 }
 Num H=0;for(int i=0;i<n;i++)H+=d[F][i];return H;
}
map<int,vector<Tile>> bank;
void collect(int target,int maxn){
 map<Key,V> prev{{0,{0}}};
 for(int n=1;n<=maxn;++n){
  map<Key,V> next;
  for(const auto& [_,r]:prev){
   Tile t=tile(r);int g=t.h;
   if(strong(r)&&g>1 && (g==9||g==27||g==81))bank[g].push_back(t);
   if(n==maxn)continue;
   int F=(1<<n)-1;
   vector<Num> cut(F+1);
   Num out[12]={};
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) out[i]+=t.through[i][j].e;
   for(unsigned a=0;a<=F;++a){
    if(a){
      int i=countr_zero(a);unsigned p=a&(a-1);
      cut[a]=cut[p]+out[i];
      for(unsigned b=p;b;b&=b-1){
        int j=countr_zero(b);
        cut[a]-=t.through[i][j].e+t.through[j][i].e;
      }
    }
    if(g+2*cut[a]>target) continue;
    V q=r;for(int i=0;i<n;++i)if(!(a>>i&1))q[i]|=1u<<n;
    q.push_back(a);next.emplace(canon(q),std::move(q));
   }
  }
  cerr<<"types "<<n<<" "<<prev.size()<<endl;
  prev.swap(next);
 }
}

constexpr int BUF=40;
int maxslots;
long long layouts=0, xhits=0, tips=0;
void inspect(const vector<const Tile*>& a){
 int l=a.size(),n=0,pos[BUF]={};
 Num pre[BUF]={1},suf[BUF]={},cover[BUF]={};suf[l]=cover[l]=1;
 for(int i=0;i<l;i++){pre[i+1]=pre[i]*a[i]->h;n+=a[i]->r.size();pos[i+1]=n;}
 for(int i=l-1;i>=0;i--){suf[i]=suf[i+1]*a[i]->h;cover[i]=cover[i+1]*a[i]->R;}
 if(pre[l]!=729||n>=32)abort();
 if(cover[0]<12217)return;
 layouts++;
 Num W[BUF][BUF]={}, starts[BUF]={}, rem[BUF]={};
 for(int i=0;i<l;i++){
  Num tail=cover[i+1];
  for(int u=pos[i];u<pos[i+1];u++){
   P st=a[i]->start[u-pos[i]];
   starts[u]=(st.e+st.o)*tail*pre[i];
   for(int v=pos[i];v<pos[i+1];v++)W[u][v]=a[i]->through[u-pos[i]][v-pos[i]].e*pre[i]*suf[i+1];
  }
  P mid{1,0};
  for(int j=i+1;j<l;j++){
   for(int u=pos[i];u<pos[i+1];u++)for(int v=pos[j];v<pos[j+1];v++)
    W[u][v]=(a[i]->start[u-pos[i]]*mid*a[j]->end[v-pos[j]]).e*pre[i]*suf[j+1];
   mid=mid*P{a[j]->R/2,a[j]->R/2};
  }
 }
 for(int j=n-1;j>=0;j--)rem[j]=rem[j+1]+starts[j];
 vector<Num> xs;for(Num x=12217;x<=rem[0]+729;x*=2)xs.push_back(x);
 auto walk=[&](auto&&self,int i,unsigned A,unsigned B,Num X,Num c)->void{
  if(c>2645)return;
  auto it=lower_bound(xs.begin(),xs.end(),X);
  if(it==xs.end()||*it>X+rem[i])return;
  if(i==pos[1]&&!A)return;
  if(i==n){
    tips++;if(*it==X)xhits++;
    if(c!=2645 || *it!=X || !(B>>pos[l-1]))return;
    V r(n+1);unsigned F=(1u<<n)-1;
    for(int j=0;j<l;j++){
     unsigned later=F^((1u<<pos[j+1])-1);
     for(int u=pos[j];u<pos[j+1];u++)r[u]=(a[j]->r[u-pos[j]]<<pos[j])|later|(((B>>u)&1)<<n);
    }
    r[n]=A;cout<<"FOUND X "<<X<<" a";for(unsigned z:r)cout<<" "<<z;cout<<endl;exit(0);
  }
  Num Aq=0,Bq=0;
  for(unsigned b=A;b;b&=b-1)Aq+=W[countr_zero(b)][i];
  for(unsigned b=B;b;b&=b-1)Bq+=W[i][countr_zero(b)];
  self(self,i+1,A|(1u<<i),B,X+starts[i],c+Bq);
  self(self,i+1,A,B|(1u<<i),X,c+Aq);
 };
 walk(walk,0,0,0,729,0);
}
void layout(const vector<const Tile*>& parts,int i,vector<const Tile*>& a,int remaining,const Tile& one){
 if(i==int(parts.size())){
  size_t p=a.size();
  for(int j=0;j<=remaining;j++){inspect(a);a.push_back(&one);}
  a.resize(p);return;
 }
 size_t p=a.size();
 for(int j=0;j<=remaining;j++){a.push_back(parts[i]);layout(parts,i+1,a,remaining-j,one);a.pop_back();a.push_back(&one);}
 a.resize(p);
}
int main(int argc, char**argv){
 int maxn=atoi(argv[1]);maxslots=atoi(argv[2]);
 collect(81,maxn);

 Tile one=tile({0}); vector<const Tile*> v;
 cout<<"BANK";

 for(int h:{9,27,81})cout<<" "<<h<<":"<<bank[h].size();
 cout<<endl;
 for(vector<int> parts: {vector<int>{9,81}, {81,9}, {27,27}, {9,9,9}}){
  vector<const Tile*> a;
  auto gen=[&](auto&& self, int i)->void{
   if(i==int(parts.size())){layout(a,0,v,maxslots,one);return;}
   for(auto&t:bank[parts[i]]){a.push_back(&t);self(self,i+1);a.pop_back();}
  };
  gen(gen,0);
 }
 cout<<"DONE "<<layouts<<" tips "<<tips<<" xhits "<<xhits<<endl;
}
