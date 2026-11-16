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
 array<V,10> groups;
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
 int n=a.size(),F=(1<<n)-1,d[1024][10];
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
   int g=quick(r);
   if(strong(r)&&g>1)bank[g].push_back(tile(r));
   if(n==maxn)continue;
   for(unsigned a=0;a<(1u<<n);++a){
    V q=r;for(int i=0;i<n;++i)if(!(a>>i&1))q[i]|=1u<<n;
    q.push_back(a);if(quick(q)<=target)next.emplace(canon(q),std::move(q));
   }
  }
  cerr<<"types "<<n<<" "<<prev.size()<<endl;
  prev.swap(next);
 }
}
constexpr int BUF=32;
int bound,limit,normalized_odd;
map<pair<int,int>,pair<Num,V>> rows;
vector<int> offsets;
long long layouts_checked=0;
void inspect(const vector<const Tile*>& a,int g){
 int l=a.size(),n=0,pos[BUF]={};Num pre[BUF]={1},suf[BUF]={};suf[l]=1;
 for(int i=0;i<l;i++){pre[i+1]=pre[i]*a[i]->h;n+=a[i]->r.size();pos[i+1]=n;}
 for(int i=l-1;i>=0;i--)suf[i]=suf[i+1]*a[i]->h;
 if(pre[l]!=g||n>=BUF)abort();
 layouts_checked++;
 Num W[BUF][BUF]={}, starts[BUF]={}, rem[BUF]={};
 for(int i=0;i<l;i++){
  Num tail=1;for(int j=i+1;j<l;j++)tail*=a[j]->R;
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
 vector<Num> slopes;for(Num k=normalized_odd;k<=rem[0];k*=2)slopes.push_back(k);
 int lowK=(g-offsets.back())/2,maxK=(g-offsets.front())/2;
 auto record=[&](unsigned A,Num S,Num K){
   if(K<lowK||!binary_search(slopes.begin(),slopes.end(),S))return;
   int d=g-2*K;if(!binary_search(offsets.begin(),offsets.end(),d))return;
   auto key=make_pair(g,d);if(rows.count(key)&&rows[key].first<=S)return;
   unsigned F=(1u<<n)-1;
   V r(n+1);
   for(int i=0;i<l;i++){
    unsigned later=(F^((1u<<pos[i+1])-1));
    for(int u=pos[i];u<pos[i+1];u++)
     r[u]=(a[i]->r[u-pos[i]]<<pos[i])|later|((!(A>>u&1))<<n);
   }
   r[n]=A;rows[key]={S,r};
 };
 auto walk=[&](auto&&self,int i,unsigned A,Num S,Num K)->void{
  if(K>maxK)return;
  if(i==pos[1]&&!A)return;
  auto t=lower_bound(slopes.begin(),slopes.end(),S);
  if(t==slopes.end()||*t>S+rem[i])return;
  if(i==n){record(A,S,K);return;}
  Num add=0;for(unsigned b=A;b;b&=b-1){int u=countr_zero(b);add+=W[u][i]+W[i][u];}
  self(self,i+1,A|(1u<<i),S+starts[i],K+add);
  if(i!=n-1||a.back()->r.size()!=1)self(self,i+1,A,S,K);
 };
 walk(walk,0,0,0,0);
}
void layout(const vector<const Tile*>& parts,int i,vector<const Tile*>& a,int g,int remaining,const Tile* one){
 if(i==int(parts.size())){
  // No trailing singleton gap needed: initially it may contain the last outcontact
  // and hence is included below, with the final bit forced to be an outcontact.
  size_t p=a.size();
  for(int j=0;j<=remaining;j++){inspect(a,g);a.push_back(one);}
  a.resize(p);return;
 }
 size_t p=a.size();
 for(int j=0;j<=remaining;j++){a.push_back(parts[i]);layout(parts,i+1,a,g,remaining-j,one);a.pop_back();a.push_back(one);}
 a.resize(p);
}
void factor(vector<const Tile*>& parts,int g,const Tile* one){
 if(g>=limit){
  int K=(g-offsets.front())/2, slots=3+(63-countl_zero((unsigned long long)K));
  vector<const Tile*> a;layout(parts,0,a,g,slots,one);
 }
 for(auto& [h,ts]:bank)if(h*g<=bound)for(auto& t:ts){
   parts.push_back(&t);factor(parts,h*g,one);parts.pop_back();
 }
}
int main(int argc,char**argv){
 if(argc!=4)return 2;
 limit=atoi(argv[1]);bound=atoi(argv[2]);normalized_odd=atoi(argv[3]);
 if(limit<3||bound>100||bound<limit||normalized_odd<1)return 2;
 offsets={-45,-43,-41,-39,-37,-35,-33,-31,-29};
 collect(bound,bound>=75?9:8);
 auto one=tile({0});vector<const Tile*> v;factor(v,1,&one);
 cout<<"LAYOUTS "<<layouts_checked<<endl;
 for(const auto&[k,t]:rows){
  auto[g,d]=k;
  cout<<g<<" "<<d<<" "<<t.first;for(unsigned r:t.second)cout<<" "<<r;cout<<"\n";
 }
}
