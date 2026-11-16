#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>
using namespace std;
using V=vector<int>;using ll=long long;
using P=array<ll,2>;
P mul(P a,P b){return {a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]};}
P add(P a,P b){return {a[0]+b[0],a[1]+b[1]};}
P sc(P a,ll b){return {a[0]*b,a[1]*b};}
struct Tile{
 V a;int n;ll h;
 P opt;
 vector<P> first,last;vector<vector<P>> full;
 Tile(V aa):a(aa),n(aa.size()),first(n),last(n),full(n,vector<P>(n)){
  int N=1<<n,F=N-1;
  vector<V> D(N,V(n));V H(N);H[0]=1;
  for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){
   int t=S^(1<<i),k=!t;for(int j=0;j<n;j++)if((t&a[i])>>j&1)k+=D[t][j];
   D[S][i]=k;H[S]+=k;
  }
  h=H[F];
  for(int i=0;i<n;i++){
   vector<V> Q(N,V(n));
   for(int S=1;S<N;S++)if(S>>i&1)for(int j=0;j<n;j++)if(S>>j&1){
    int v=0;
    if(i==j)v=(S==(1<<i));
    else{int t=S^(1<<j);for(int k=0;k<n;k++)if((t&~a[j])>>k&1)v+=Q[t][k];}
    Q[S][j]=v;ll w=ll(v)*H[F^S];int p=__builtin_popcount(S)%2;
    first[i][p]+=w;last[j][p]+=w;full[i][j][p]+=w;
   }
  }
  opt={h,0};for(P z:first)opt=add(opt,z);
 }
};
bool strong(const V&a){int n=a.size(),F=(1<<n)-1;for(int rev=0;rev<2;rev++){
 int m=1,q=1;while(q){int i=__builtin_ctz(q);q&=q-1;int b=0;if(!rev)b=a[i];else for(int j=0;j<n;j++)if(a[j]>>i&1)b|=1<<j;q|=b&~m;m|=b;}if(m!=F)return false;}return true;}
uint32_t canon(const V&a){V p(a.size());for(int i=0;i<(int)p.size();i++)p[i]=i;uint32_t best=UINT32_MAX;do{uint32_t r=0;for(int i=0;i<(int)p.size();i++)for(int j=i+1;j<(int)p.size();j++)r=2*r+(a[p[i]]>>p[j]&1);best=min(best,r);}while(next_permutation(p.begin(),p.end()));return best;}
map<int,map<pair<int,uint32_t>,Tile>> bank;
int dp[256][8];
int total(const V&a){
 int n=a.size(),F=(1<<n)-1; int h=0;
 for(int s=1;s<=F;s++){h=0;for(int b=s;b;b&=b-1){
  int i=__builtin_ctz(b),r=s^(1<<i),v=!r;
  for(int x=r&a[i];x;x&=x-1)v+=dp[r][__builtin_ctz(x)];
  dp[s][i]=v;h+=v;
 }}return h;
}
int max_tile,goal; int divisor(int q){
 return q>0&&goal%q==0;
}
void gen(V a){
 int h=total(a);if(h>goal)return;
 if(h>1&&strong(a))bank[h].emplace(make_pair(a.size(),canon(a)),Tile(a));
 int n=a.size();if(n==max_tile)return;
 for(int b=0;b<(1<<(n-1));b++){V c=a;for(int i=0;i<n;i++)if(!(b>>i&1))c[i]|=1<<n;c.push_back(b);gen(c);}
}

int mult,S,Nmax;ll low;
ll hosts=0;
map<pair<ll,ll>,int> support;
map<int,set<ll>> cover;map<ll,int> signed_zero;
map<array<ll,3>,V> witnesses;
void inspect(const vector<const Tile*>& a){
 ++hosts;
 int k=a.size(),n=0; V off;ll g=1;
 for(auto t:a){off.push_back(n);n+=t->n;g*=t->h;}
 int F=(1<<n)-1,first=(1<<a.front()->n)-1,tail=F^((1<<off.back())-1);
 ll E[24]={},O[24]={},D[24][24]={},W[24][24]={};
 ll pref[25]={1},suf[25]={};
 suf[k]=1;
 for(int i=0;i<k;i++)pref[i+1]=pref[i]*a[i]->h;
 for(int i=k-1;i>=0;i--)suf[i]=suf[i+1]*a[i]->h;
 for(int i=0;i<k;i++){
  auto t=a[i];P after={1,0};for(int j=i+1;j<k;j++)after=mul(after,a[j]->opt);
  for(int u=0;u<t->n;u++){
   int U=off[i]+u;auto io=sc(mul(t->first[u],after),pref[i]);E[U]=io[1];O[U]=io[0];
   for(int v=0;v<t->n;v++)D[U][off[i]+v]=pref[i]*suf[i+1]*t->full[u][v][0];
   P mid={1,0};for(int j=i+1;j<k;j++){
    auto z=a[j];for(int v=0;v<z->n;v++)D[U][off[j]+v]=pref[i]*suf[j+1]*mul(mul(t->first[u],mid),z->last[v])[0];
    mid=mul(mid,z->opt);
   }
  }
 }
 ll orient[24]={};for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){orient[i]+=D[i][j];W[i][j]=D[i][j]+D[j][i];}
 ll cut=0,e=0,o=g;int prev=0;
 for(int z=1;z<=F;z++){
  int mask=z^(z>>1),toggle=prev^mask,i=__builtin_ctz(toggle);
  ll diff=orient[i];for(int t=prev&~toggle;t;t&=t-1)diff-=W[i][__builtin_ctz(t)];
  int sign=((mask&toggle)!=0)?1:-1;
  cut+=sign*diff;e+=sign*E[i];o+=sign*O[i];prev=mask;
  ll h=mult*(g+2*cut);
  if(h<low||h>low+4||!(mask&first)||!((F^mask)&tail))continue;
  ll ee=mult*e,oo=mult*o;
  support[{ee,oo}]|=1<<((h-low)/2);
  for(auto [b,RR]:cover)if(mult%b==0)for(ll R:RR){
   ll u=ll(mult/b)*R*(e+o);
   signed_zero[u/(u&-u)]|=1<<((h-low)/2);
  }
  auto key=array<ll,3>{ee,oo,h};
  if(!witnesses.count(key)){
   V r(n+1);for(int l=0;l<k;l++)for(int u=0;u<a[l]->n;u++){
    int U=off[l]+u;r[U]=(a[l]->a[u]<<off[l])|(F^((1<<(off[l]+a[l]->n))-1));
    if(!(mask>>U&1))r[U]|=1<<n;
   }r[n]=mask;r.push_back(mult);
   witnesses.emplace(key,r);
  }
 }
}
Tile one({0});
void layouts(const vector<const Tile*>& vec,int i,vector<const Tile*> out,int n){
 if(i==(int)vec.size()){for(int j=n;j<=Nmax;j++){inspect(out);out.push_back(&one);}return;}
 auto t=vec[i];for(int j=n;j+t->n<=Nmax;j++){
  auto v=out;v.push_back(t);layouts(vec,i+1,v,j+t->n);out.push_back(&one);
 }
}
void factors(V f,int i,vector<const Tile*> vv){
 if(i==(int)f.size()){if(vv.empty()){for(int j=1;j<=Nmax;j++){vv.push_back(&one);inspect(vv);}}else layouts(vv,0,{},0);return;}
 for(auto& [k,t]:bank[f[i]]){vv.push_back(&t);factors(f,i+1,vv);vv.pop_back();}
}
map<int,vector<V>> facs;
void fact(int n,V f){
 if(n==1){facs[1].push_back(f);return;}
 for(auto &[p,b]:bank)if(n%p==0&&p>1&&b.size()){auto ff=f;ff.push_back(p);fact(n/p,ff);}
}
int busch(int n){int a[]={3,5,9};int h=a[(n-3)%3];for(int j=0;j<(n-3)/3;j++)h*=5;return h;}
int main(int argc,char**argv){
 S=3;goal=73;low=756ll*(1ll<<S)-33;max_tile=8;gen({0});
 for(goal=1;goal<=73;goal+=2){
 cover.clear();support.clear();witnesses.clear();signed_zero.clear();hosts=0;
 cerr<<"goal "<<goal<<"\n";
 for(int m=1;m<=goal;m++)if(goal%m==0){
  facs.clear();fact(m,{});for(V f:facs[1]){
   set<ll> v={1};for(int p:f){set<ll>w;for(auto &[k,t]:bank[p])for(auto x:v)w.insert(x*(t.opt[0]+t.opt[1]));v=w;}
   cover[m].insert(v.begin(),v.end());
  }
 }
 cerr<<"bank";for(auto &[p,b]:bank)cerr<<" "<<p<<":"<<b.size();cerr<<"\n";
 for(int m=1;m<=goal;m++)if(goal%m==0){
  mult=m;bool ok=false;for(int j=0;j<=4;j+=2)ok|=(low+j)%m==0;if(!ok)continue;
  Nmax=1;for(int n=3;n<=23;n++){if(ll(busch(n))*m>low+4)break;Nmax=n-1;}
  if(Nmax>22)return 2;
  facs.clear();fact(goal/m,{});for(V f:facs[1])factors(f,0,{});
  cerr<<"mult "<<m<<" deletion ceiling "<<Nmax<<" cumulative hosts "<<hosts<<"\n";
 }
 for(auto [odd,b]:signed_zero)if(b==7)cerr<<"signzero "<<S<<" "<<odd<<"\n";
 for(auto [eo,b]:support){auto [e,o]=eo;cout<<goal<<" "<<e<<" "<<o<<" "<<b<<"\n";}
 for(auto& [key,r]:witnesses){auto [e,o,h]=key;if(support[{e,o}]!=7)continue;cerr<<"triple "<<S<<" "<<h<<" "<<e<<" "<<o<<" mult "<<r.back()<<" rows";for(int j=0;j+1<(int)r.size();j++)cerr<<" "<<r[j];cerr<<"\n";}
}
}
