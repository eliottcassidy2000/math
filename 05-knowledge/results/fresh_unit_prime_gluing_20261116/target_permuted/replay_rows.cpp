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
void gen(V a){Tile t(a);if(t.h>15)return;if(t.h>1&&(t.h==3||t.h==5||t.h==15)&&strong(a))bank[t.h].emplace(make_pair(a.size(),canon(a)),t);
 int n=a.size();if(n==6)return;for(int b=0;b<(1<<(n-1));b++){V c=a;for(int i=0;i<n;i++)if(!(b>>i&1))c[i]|=1<<n;c.push_back(b);gen(c);}}
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
int main(int argc,char**argv){
 if(argc!=2)return 2;S=stoi(argv[1]);if(S<1||S>8)return 2;low=756ll*(1ll<<S)-33;gen({0});
 cover[1]={1};for(int g:{3,5,15})for(auto& [k,t]:bank[g])cover[g].insert(t.opt[0]+t.opt[1]);
 for(auto a:cover[3])for(auto b:cover[5])cover[15].insert(a*b);
 cerr<<"bank "<<bank[3].size()<<" "<<bank[5].size()<<" "<<bank[15].size()<<"\n";
 for(int m:{1,3,5,15}){
  mult=m;bool ok=false;for(int j=0;j<=4;j+=2)ok|=(low+j)%m==0;if(!ok)continue;
  int f[]={3,5,9};Nmax=1;for(int n=3;n<=23;n++){if(ll(f[(n-3)%3])*m>low+4)break;Nmax=n-1;if(n%3==2)for(int &t:f)t*=5;}
  if(Nmax>22)return 2;
  if(m==1){factors({15},0,{});factors({3,5},0,{});factors({5,3},0,{});}
  if(m==3)factors({5},0,{});
  if(m==5)factors({3},0,{});
  if(m==15)factors({},0,{});
  cerr<<"mult "<<m<<" deletion ceiling "<<Nmax<<" cumulative hosts "<<hosts<<"\n";
 }
 for(auto [odd,b]:signed_zero)if(b==7)cerr<<"signzero "<<S<<" "<<odd<<"\n";
 for(auto [eo,b]:support){auto [e,o]=eo;cout<<e<<" "<<o<<" "<<b<<"\n";}
 for(auto& [key,r]:witnesses){auto [e,o,h]=key;if(__builtin_popcount(support[{e,o}])<2)continue;cerr<<"triple "<<S<<" "<<h<<" "<<e<<" "<<o<<" mult "<<r.back()<<" rows";for(int j=0;j+1<(int)r.size();j++)cerr<<" "<<r[j];cerr<<"\n";}
}
